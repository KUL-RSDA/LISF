!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_qclaisoilm
! \label{NoahMP401_qclaisoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 19 Jan 2022: Samuel Scherrer: Modified from da_soilm for da_laisoilm
!
! !INTERFACE:
subroutine NoahMP401_qclaisoilm(n, LSM_State)

    ! !USES:
    use ESMF
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_surface
    use LIS_logMod,  only  : LIS_verify, LIS_logunit
    use NoahMP401_lsmMod
    use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in)    :: n
    type(ESMF_State)       :: LSM_State
    !
    ! !DESCRIPTION:
    !
    !  Returns the LAI and soilmoisture related state prognostic variables for
    !  data assimilation
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[n] index of the nest \newline
    !  \item[LSM\_State] ESMF State container for LSM state variables \newline
    !  \end{description}
    !EOP
    integer                :: t, status
    type(ESMF_Field)       :: laiField
    type(ESMF_Field)       :: sm1Field
    real, pointer          :: lai(:)
    real, pointer          :: soilm1(:)
    real                   :: laimax
    real                   :: laimin
  integer                :: gid
  real                   :: laitmp

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: laimean(LIS_rc%ngrid(n))
  integer                :: nlaimean(LIS_rc%ngrid(n))
 
  integer                :: N_ens
  real                   :: state_tmp(LIS_rc%nensem(n)),state_mean


    real                   :: smmax1!,smmax2,smmax3,smmax4
    real                   :: smmin1!,smmin2,smmin3,smmin4

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in NoahMP401_qcsoilm")
 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in NoahMP401_qcsoilm")

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in NoahMP401_qcsoilm")

  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in NoahMP401_qcsoilm")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(soilm1(t).gt.smmax1) soilm1(t) = smmax1
     if(soilm1(t).lt.smmin1) soilm1(t) = smmin1
  enddo


  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(laiField,"Max Value",laimax,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(laiField,"Min Value",laimin,rc=status)
  call LIS_verify(status)



  update_flag    = .true.
  perc_violation = 0.0
  laimean       = 0.0
  nlaimean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t)

     if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
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
           if((lai(t).gt.laimin).and.&
                (lai(t).lt.laimax)) then 
              laimean(gid) = laimean(gid) + &
                   lai(t) 
              nlaimean(gid) = nlaimean(gid) + 1
           endif
        endif
     endif
  enddo
  
  do gid=1,LIS_rc%ngrid(n)
     if(nlaimean(gid).gt.0) then
        laimean(gid) = laimean(gid)/nlaimean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     laitmp =  lai(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        lai(t) = laitmp
     elseif(perc_violation(gid).lt.0.8) then
        if(laitmp.lt.laimin.or.laitmp.gt.laimax) then
           lai(t) = laimean(gid)
        else
           lai(t) = lai(t) 
        endif
     else
         lai(t) = noahmp401_struc(n)%noahmp401(t)%lai
     endif
  enddo


end subroutine NoahMP401_qclaisoilm

