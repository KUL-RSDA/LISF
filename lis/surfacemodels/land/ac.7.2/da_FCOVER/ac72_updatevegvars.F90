!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: ac72_updatevegvars
! \label{ac72_updatevegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac72_updatevegvars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use ac72_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)       :: AC72CCiprevField, AC72CCiprevIncrField
  type(ESMF_Field)       :: AC72BiomassField, AC72BiomassIncrField

  integer                :: t,gid
  integer                :: status
  real, pointer          :: AC72CCiprev(:), AC72CCiprevincr(:)
  real, pointer          :: AC72Biomass(:), AC72Biomassincr(:)
  real                   :: AC72CCiprevtmp,AC72CCiprevmax,AC72CCiprevmin
  real                   :: AC72Biomasstmp,AC72Biomassmax,AC72Biomassmin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: AC72CCiprevmean(LIS_rc%ngrid(n))
  real                   :: AC72Biomassmean(LIS_rc%ngrid(n))
  integer                :: nAC72CCiprevmean(LIS_rc%ngrid(n))
  integer                :: nAC72Biomassmean(LIS_rc%ngrid(n))

 
  call ESMF_StateGet(LSM_State,"AC72 CCiprev",AC72CCiprevField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"AC72 CCiprev",AC72CCiprevIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(AC72CCiprevField,localDE=0,farrayPtr=AC72CCiprev,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(AC72CCiprevIncrField,localDE=0,farrayPtr=AC72CCiprevincr,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(AC72CCiprevField,"Max Value",AC72CCiprevmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(AC72CCiprevField,"Min Value",AC72CCiprevmin,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LSM_State,"AC72 Biomass",AC72BiomassField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"AC72 Biomass",AC72BiomassIncrField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(AC72BiomassField,localDE=0,farrayPtr=AC72Biomass,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(AC72BiomassIncrField,localDE=0,farrayPtr=AC72Biomassincr,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(AC72BiomassField,"Max Value",AC72Biomassmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(AC72BiomassField,"Min Value",AC72Biomassmin,rc=status)
  call LIS_verify(status)


  update_flag = .true.
  perc_violation = 0.0
  AC72CCiprevmean = 0.0
  nAC72CCiprevmean = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72CCiprevtmp =  AC72CCiprev(t) + AC72CCiprevincr(t)


     if(AC72CCiprevtmp.lt.AC72CCiprevmin.or.AC72CCiprevtmp.gt.AC72CCiprevmax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) + 1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((AC72CCiprev(t)+AC72CCiprevincr(t).gt.AC72CCiprevmin).and.&
                (AC72CCiprev(t)+AC72CCiprevincr(t).lt.AC72CCiprevmax)) then 
              AC72CCiprevmean(gid) = AC72CCiprevmean(gid) + &
                   AC72CCiprev(t) + AC72CCiprevincr(t)
              nAC72CCiprevmean(gid) = nAC72CCiprevmean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nAC72CCiprevmean(gid).gt.0) then
        AC72CCiprevmean(gid) = AC72CCiprevmean(gid)/nAC72CCiprevmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72CCiprevtmp =  AC72CCiprev(t) + AC72CCiprevincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        AC72CCiprev(t) = AC72CCiprevtmp
     elseif(perc_violation(gid).lt.0.8) then
        if(AC72CCiprevtmp.lt.AC72CCiprevmin.or.AC72CCiprevtmp.gt.AC72CCiprevmax) then
           AC72CCiprev(t) = AC72CCiprevmean(gid)
        else
           AC72CCiprev(t) = AC72CCiprev(t) + AC72CCiprevincr(t)
        endif
     endif
  enddo


  update_flag = .true.
  perc_violation = 0.0
  AC72Biomassmean = 0.0
  nAC72Biomassmean = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72Biomasstmp =  AC72Biomass(t) + AC72Biomassincr(t)


     if(AC72Biomasstmp.lt.AC72Biomassmin.or.AC72Biomasstmp.gt.AC72Biomassmax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) + 1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((AC72Biomass(t)+AC72Biomassincr(t).gt.AC72Biomassmin).and.&
                (AC72Biomass(t)+AC72Biomassincr(t).lt.AC72Biomassmax)) then 
              AC72Biomassmean(gid) = AC72Biomassmean(gid) + &
                   AC72Biomass(t) + AC72Biomassincr(t)
              nAC72Biomassmean(gid) = nAC72Biomassmean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nAC72Biomassmean(gid).gt.0) then
        AC72Biomassmean(gid) = AC72Biomassmean(gid)/nAC72Biomassmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72Biomasstmp =  AC72Biomass(t) + AC72Biomassincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        AC72Biomass(t) = AC72Biomasstmp
     elseif(perc_violation(gid).lt.0.8) then
        if(AC72Biomasstmp.lt.AC72Biomassmin.or.AC72Biomasstmp.gt.AC72Biomassmax) then
           AC72Biomass(t) = AC72Biomassmean(gid)
        else
           AC72Biomass(t) = AC72Biomass(t) + AC72Biomassincr(t)
        endif
     endif
  enddo
end subroutine ac72_updatevegvars