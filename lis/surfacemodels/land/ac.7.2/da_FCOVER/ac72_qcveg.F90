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
! !ROUTINE: ac72_qcveg
! \label{ac72_qcveg}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
! !INTERFACE:
subroutine ac72_qcveg(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use ac72_lsmMod
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  QC's the related state prognostic variable objects for
!  veg data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: AC72CCiprevField
  integer                :: t
  integer                :: status
  real, pointer          :: AC72CCiprev(:)

  real                   :: AC72CCiprevmax
  real                   :: AC72CCiprevmin

  integer                :: gid
  real                   :: AC72CCiprevtmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: AC72CCiprevmean(LIS_rc%ngrid(n))
  integer                :: nAC72CCiprevmean(LIS_rc%ngrid(n))
 
  integer                :: N_ens
  real                   :: state_tmp(LIS_rc%nensem(n)),state_mean

  call ESMF_StateGet(LSM_State,"AC72 CCiprev",AC72CCiprevField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(AC72CCiprevField,localDE=0,farrayPtr=AC72CCiprev,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(AC72CCiprevField,"Max Value",AC72CCiprevmax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(AC72CCiprevField,"Min Value",AC72CCiprevmin,rc=status)
  call LIS_verify(status)


  update_flag    = .true.
  perc_violation = 0.0
  AC72CCiprevmean = 0.0
  nAC72CCiprevmean = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72CCiprevtmp =  AC72CCiprev(t)

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
           if((AC72CCiprev(t).gt.AC72CCiprevmin).and.&
                (AC72CCiprev(t).lt.AC72CCiprevmax)) then 
              AC72CCiprevmean(gid) = AC72CCiprevmean(gid) + &
                   AC72CCiprev(t) 
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

     AC72CCiprevtmp =  AC72CCiprev(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        AC72CCiprev(t) = AC72CCiprevtmp
     elseif(perc_violation(gid).lt.0.8) then
        if(AC72CCiprevtmp.lt.AC72CCiprevmin.or.AC72CCiprevtmp.gt.AC72CCiprevmax) then
           AC72CCiprev(t) = AC72CCiprevmean(gid)
        else
           AC72CCiprev(t) = AC72CCiprev(t) 
        endif
     endif
  enddo

#if 0 
  N_ens = LIS_rc%nensem(n)
  do t=1,N_ens
     state_tmp(t) = AC72CCiprev(t)
  enddo
  state_mean =sum(state_tmp)/N_ens

  write(113,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,21F8.3)') &
       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr, &
       LIS_rc%mn, LIS_rc%ss, &
       state_mean, &
       state_tmp

#endif
end subroutine ac72_qcveg