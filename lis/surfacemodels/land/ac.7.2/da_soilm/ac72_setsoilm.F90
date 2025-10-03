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
! !ROUTINE: ac72_setsoilm
!  \label{ac72_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari: Modified for ac72 
! 18Apr 2018: Mahdi Navari: Bug fixed
! 
! Apply the update if it met the update conditions
! Update conditions: 
!                  1- Prior SM(sh2o) + increment > MIN_THRESHOLD 
!                  2- Prior SM(sh2o) + increment < sm_threshold
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states    


! !INTERFACE:
subroutine ac72_setsoilm(n, LSM_State)
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac72_lsmMod
  use ac_kinds, only: sp

  implicit none

  integer, intent(in) :: n
  type(ESMF_State)    :: LSM_State

  real, parameter     :: MIN_THRESHOLD = 0.02
  real                :: MAX_THRESHOLD, sm_threshold, tmpval
  integer             :: t, j, i, gid, m, t_unpert, status
  integer             :: RESULT, pcount, icount, nIter
  logical             :: bounds_violation
  logical             :: update_flag(LIS_rc%ngrid(n))
  logical             :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical             :: flag_ens(LIS_rc%ngrid(n))
  logical             :: flag_tmp(LIS_rc%nensem(n))
  logical             :: update_flag_ens(LIS_rc%ngrid(n))
  logical             :: update_flag_new(LIS_rc%ngrid(n))

  integer             :: SOILTYP
  integer, parameter         :: n_layers = 10  ! Replace with variable if needed elsewhere
  type(ESMF_Field), allocatable :: smField(:)
  real, pointer, dimension(:), allocatable :: soilm(:,:)

  real                :: delta(:)
  real, allocatable   :: tmp(:,:), MaxEnsSM(:), MinEnsSM(:)
  real                :: smc_tmp

  allocate(smField(n_layers))
  allocate(soilm(n_layers, :))
  allocate(tmp(n_layers, LIS_rc%nensem(n)))
  allocate(MaxEnsSM(n_layers), MinEnsSM(n_layers))
  allocate(delta(n_layers))

  ! Get soil moisture fields
  do j = 1, n_layers
     call ESMF_StateGet(LSM_State, "Soil Moisture Layer "//trim(adjustl(itoa(j))), smField(j), rc=status)
     call LIS_verify(status, "ESMF_StateGet failed for Soil Moisture Layer "//trim(adjustl(itoa(j)))//" in ac72_setsoilm")
     call ESMF_FieldGet(smField(j), localDE=0, farrayPtr=soilm(j,:), rc=status)
     call LIS_verify(status, "ESMF_FieldGet failed for Soil Moisture Layer "//trim(adjustl(itoa(j)))//" in ac72_setsoilm")
  enddo

  update_flag = .true.
  update_flag_tile = .true.

  ! Loop over tiles
  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     SOILTYP = AC72_struc(n)%ac72(t)%soiltype
     MAX_THRESHOLD = AC72_struc(n)%ac72(t)%soillayer(1)%sat/100.0 - epsilon(0._sp)
     sm_threshold = MAX_THRESHOLD - 0.02

     gid = LIS_domain(n)%gindex( &
           LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col, &
           LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row)

     do j = 1, n_layers
        delta(j) = soilm(j, t) - AC72_struc(n)%ac72(t)%smc(j)
        if (AC72_struc(n)%ac72(t)%smc(j) + delta(j) > MIN_THRESHOLD .and. &
            AC72_struc(n)%ac72(t)%smc(j) + delta(j) < sm_threshold) then
           update_flag(gid) = update_flag(gid) .and. .true.
           update_flag_tile(t) = update_flag_tile(t) .and. .true.
        else
           update_flag(gid) = update_flag(gid) .and. .false.
           update_flag_tile(t) = update_flag_tile(t) .and. .false.
        endif
     enddo
  enddo

  ! Ensemble majority filter
  update_flag_ens = .true.
  do i = 1, LIS_rc%npatch(n,LIS_rc%lsm_index), LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex(&
           LIS_surface(n, LIS_rc%lsm_index)%tile(i)%col, &
           LIS_surface(n, LIS_rc%lsm_index)%tile(i)%row)
     flag_tmp = update_flag_tile(i:i + LIS_rc%nensem(n) - 1)
     pcount = COUNT(flag_tmp)
     if (pcount < LIS_rc%nensem(n) * 0.5) then
        update_flag_ens(gid) = .false.
     endif
     update_flag_new(gid) = update_flag(gid) .or. update_flag_ens(gid)
  enddo

  ! Apply updates
  do i = 1, LIS_rc%npatch(n,LIS_rc%lsm_index), LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex( &
           LIS_surface(n, LIS_rc%lsm_index)%tile(i)%col, &
           LIS_surface(n, LIS_rc%lsm_index)%tile(i)%row)

     if (update_flag_new(gid)) then
        tmp = LIS_rc%udef
        do m = 1, LIS_rc%nensem(n)
           t = i + m - 1
           if (update_flag_tile(t)) then
              do j = 1, n_layers
                 tmp(j, m) = soilm(j, t)
              enddo
           endif
        enddo

        MaxEnsSM = -1.0e6
        MinEnsSM = 1.0e6
        do j = 1, n_layers
           do m = 1, LIS_rc%nensem(n)
              if (tmp(j, m) /= LIS_rc%udef) then
                 MaxEnsSM(j) = max(MaxEnsSM(j), tmp(j, m))
                 MinEnsSM(j) = min(MinEnsSM(j), tmp(j, m))
              endif
           enddo
        enddo

        do m = 1, LIS_rc%nensem(n)
           t = i + m - 1
           if (update_flag_tile(t)) then
              do j = 1, n_layers
                 AC72_struc(n)%ac72(t)%smc(j) = soilm(j, t)
                 if (soilm(j, t) < 0.0) then
                    print*, "Negative soilm detected at t=", t, " layer=", j, " val=", soilm(j, t)
                    stop
                 endif
              enddo
           else
              do j = 1, n_layers
                 smc_tmp = 0.5 * (MaxEnsSM(j) + MinEnsSM(j))
                 AC72_struc(n)%ac72(t)%smc(j) = smc_tmp
              enddo
           endif
        enddo
     endif
  enddo

  deallocate(smField, soilm, tmp, MaxEnsSM, MinEnsSM, delta)

end subroutine ac72_setsoilm

