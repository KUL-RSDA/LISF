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
! !ROUTINE: ac70_setsoilmLAI
!  \label{ac70_setsoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari: Modified for ac70 
! 18Apr 2018: Mahdi Navari: Bug fixed
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
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
subroutine ac70_setsoilmLAI(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac70_lsmMod
  use module_sf_noahaclsm_36  !MN
  use AC_VEG_PARAMETERS_70
  use ac_kinds, only: dp

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture and LAI prognostic variables to noah's
!  model space. 
! 
!EOP
  real, parameter        :: MIN_THRESHOLD = 0.02 
  real                   :: MAX_threshold
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: AC70BIOMASSField
  real, pointer          :: soilm1(:)
  real, pointer          :: AC70BIOMASS(:)
  integer                :: t, j,i, gid, m, t_unpert
  integer                :: status
  real                   :: delta(1)
  real                   :: delta1
  real                   :: tmpval
  logical                :: bounds_violation
  integer                :: nIter
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
! mn
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: tmp(LIS_rc%nensem(n)), tmp0(LIS_rc%nensem(n))
  real                   :: tmp1(LIS_rc%nensem(n))
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: flag_ens(LIS_rc%ngrid(n))
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  integer                :: RESULT, pcount, icount
  real                   :: MaxEnsSM1 
  real                   :: MinEnsSM1 
  real                   :: smc_rnd, smc_tmp 
  INTEGER, DIMENSION (1) :: seed 

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in ac70_setsoilmLAI")
  call ESMF_StateGet(LSM_State,"AC70 BIOMASS",AC70BIOMASSField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: AC70BIOMASS failed in ac70_setsoilmLAI")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac70_setsoilmLAI")
  call ESMF_FieldGet(AC70BIOMASSField,localDE=0,farrayPtr=AC70BIOMASS,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC70BIOMASS failed in ac70_setsoilmLAI")

  update_flag = .true. 
  update_flag = .true. 
  update_flag_tile= .true. 
  update_flag_tile= .true. 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  
 ! MN: NOTE: SMCMAX and SMCWLT are not stored in the data structure but we 
 !       can get module variables MAXSMC and WLTSMC from the module_sf_noahaclsm_36       
     SOILTYP = AC70_struc(n)%ac70(t)%soiltype        
     MAX_THRESHOLD = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - epsilon(0._dp)
     sm_threshold = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - 0.02
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     !MN: delta = X(+) - X(-)
     !NOTE: "ac70_updatesoilm.F90" updates the soilm_(t)     
     delta1 = soilm1(t)-AC70_struc(n)%ac70(t)%smc(1)

     ! MB: AC70
     ! MN: check    MIN_THRESHOLD < volumetric liquid soil moisture < threshold 
     if(AC70_struc(n)%ac70(t)%smc(1)+delta1.gt.MIN_THRESHOLD .and.&
          AC70_struc(n)%ac70(t)%smc(1)+delta1.lt.&
          sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   (64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
  enddo

  !!! MB: AC70
  update_flag_ens = .True.
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
      flag_tmp=update_flag_tile(i:i+LIS_rc%nensem(n)-1)
      !flag_tmp=update_flag_tile((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n))
      pcount = COUNT(flag_tmp) ! Counts the number of .TRUE. elements
      if (pcount.lt.LIS_rc%nensem(n)*0.5) then   ! 50%
         update_flag_ens(gid)= .False.
      endif
      update_flag_new(gid)= update_flag(gid).or.update_flag_ens(gid)  ! new flag
  enddo 
  
  ! MN print 
  if(i.eq.66) then !i=66  ! --> domain's center  1376
  if(LIS_rc%hr.eq.12) then
     write(2001,'(I4, 2x, 3(I2,x), 2x, 23(L1,2x))'),&
          i, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,update_flag_tile&
          ((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n)),&
          update_flag_ens(i), update_flag_new(i), update_flag(i) 
  endif !mn
  endif
  
!!! MB: AC70

  ! update step
  ! loop over grid points 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
     
     !if(update_flag(gid)) then
     if(update_flag_new(gid)) then 
!-----------------------------------------------------------------------------------------
       ! 1- update the states
       ! 1-1- if update flag for tile is TRUE --> apply the DA update    
       ! 1-2- if update flag for tile is FALSE --> resample the states  
       ! 2- adjust the states
!-----------------------------------------------------------------------------------------
        ! store update value for  cases that update_flag_tile & update_flag_new are TRUE
        ! update_flag_tile = TRUE --> means met the both min and max threshold
      
        tmp1 = LIS_rc%udef
	!icount = 1
        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           !t = (i-1)*LIS_rc%nensem(n)+m
 
	  if(update_flag_tile(t)) then
          
           tmp1(m) = soilm1(t) 
           !icount = icount + 1 
          endif
        enddo
        
        MaxEnsSM1 = -10000

        MinEnsSM1 = 10000

        do m=1,LIS_rc%nensem(n)
           if(tmp1(m).ne.LIS_rc%udef) then 
              MaxEnsSM1 = max(MaxEnsSM1, tmp1(m))

              MinEnsSM1 = min(MinEnsSM1, tmp1(m))
              
           endif
        enddo

        
        ! loop over tile       
        do m=1,LIS_rc%nensem(n)
           t = i+m-1           
           !t = (i-1)*LIS_rc%nensem(n)+m
           
           ! MN check update status for each tile  
           if(update_flag_tile(t)) then
              
              delta1 = soilm1(t)-AC70_struc(n)%ac70(t)%smc(1)

              AC70_struc(n)%ac70(t)%smc(1) = soilm1(t)
              if(soilm1(t).lt.0) then 
                 print*, 'ac70setsoilm1 ',t,soilm1(t)
                 stop
              endif
!-----------------------------------------------------------------------------------------              
              ! randomly resample the smc from [MIN_THRESHOLD,  Max value from DA @ that tiem step]
!-----------------------------------------------------------------------------------------
           else 
          
!-----------------------------------------------------------------------------------------  
! set the soil moisture to the ensemble mean  
!-----------------------------------------------------------------------------------------
              
              ! Assume sh2o = smc (i.e. ice content=0) 
              smc_tmp = (MaxEnsSM1 - MinEnsSM1)/2 + MinEnsSM1
              AC70_struc(n)%ac70(t)%smc(1) = smc_tmp
                          
           endif ! flag for each tile

        enddo ! loop over tile
       
     else ! if update_flag_new(gid) is FALSE   
        if(LIS_rc%pert_bias_corr.eq.1) then           
           !--------------------------------------------------------------------------
           ! if no update is made, then we need to readjust the ensemble if pert bias
           ! correction is turned on because the forcing perturbations may cause 
           ! biases to persist. 
           !--------------------------------------------------------------------------
           bounds_violation = .true. 
           nIter = 0
           ens_flag = .true. 
           
           do while(bounds_violation) 
              niter = niter + 1
              !t_unpert = i*LIS_rc%nensem(n)
	          t_unpert = i+LIS_rc%nensem(n)-1
              do j=1,1
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (AC70_struc(n)%ac70(t)%smc(j) - &
                            AC70_struc(n)%ac70(t_unpert)%smc(j))
                    endif
                    
                 enddo
              enddo
              
              do j=1,1
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    MAX_THRESHOLD = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - epsilon(0._dp)
                    sm_threshold = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - 0.02
                    
                    tmpval = AC70_struc(n)%ac70(t)%smc(j) - &
                         delta(j)
                    if(tmpval.le.MIN_THRESHOLD) then 
                       AC70_struc(n)%ac70(t)%smc(j) = &
                            max(AC70_struc(n)%ac70(t_unpert)%smc(j),&
                            MIN_THRESHOLD)
                       ens_flag(m) = .false. 
                    elseif(tmpval.ge.sm_threshold) then
                       AC70_struc(n)%ac70(t)%smc(j) = &
                            min(AC70_struc(n)%ac70(t_unpert)%smc(j),&
                            sm_threshold)
                       ens_flag(m) = .false. 
                    endif
                 enddo
              enddo
              
              !--------------------------------------------------------------------------
              ! Recalculate the deltas and adjust the ensemble
              !--------------------------------------------------------------------------
              do j=1,1
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (AC70_struc(n)%ac70(t)%smc(j) - &
                            AC70_struc(n)%ac70(t_unpert)%smc(j))
                    endif
                 enddo
              enddo
              
              do j=1,1
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(ens_flag(m)) then 
                       tmpval = AC70_struc(n)%ac70(t)%smc(j) - &
                            delta(j)
                       MAX_THRESHOLD = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - epsilon(0._dp)
                       sm_threshold = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - 0.02
                       
                       if(.not.(tmpval.le.0.0 .or.&
                            tmpval.gt.(MAX_THRESHOLD))) then 
                          
                          AC70_struc(n)%ac70(t)%smc(j) = &
                               AC70_struc(n)%ac70(t)%smc(j) - delta(j)
                          bounds_violation = .false.
                       endif
                    endif
                    
                    tmpval = AC70_struc(n)%ac70(t)%smc(j)
                    MAX_THRESHOLD = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - epsilon(0._dp)
                    sm_threshold = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - 0.02
                    
                    if(tmpval.le.0.0 .or.&
                         tmpval.gt.(MAX_THRESHOLD)) then 
                       bounds_violation = .true. 
                    else
                       bounds_violation = .false.
                    endif
                 enddo
              enddo
              
              if(nIter.gt.10.and.bounds_violation) then 
                 !--------------------------------------------------------------------------
                 ! All else fails, set to the bounds
                 !--------------------------------------------------------------------------
                 
                 write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
                 write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'

                 do j=1,1
                    do m=1,LIS_rc%nensem(n)
                       t = i+m-1
                       !t = (i-1)*LIS_rc%nensem(n)+m
                       
                       MAX_THRESHOLD = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - epsilon(0._dp)
                       sm_threshold = AC70_struc(n)%ac70(t)%soillayer(1)%sat/100.0 - 0.02
                       
                       if(AC70_struc(n)%ac70(t)%smc(j).gt.MAX_THRESHOLD) then
                          AC70_struc(n)%ac70(t)%smc(j) = MAX_THRESHOLD
                       endif
                       
                       if(AC70_struc(n)%ac70(t)%smc(j).lt.MIN_THRESHOLD) then 
                          AC70_struc(n)%ac70(t)%smc(j) = MIN_THRESHOLD
                       endif
!                    print*, i, m
!                    print*, 'smc',t, AC70_struc(n)%ac70(t)%smc(:)
!                    print*, 'sh2o ',t,AC70_struc(n)%ac70(t)%sh2o(:)
!                    print*, 'max ',t,MAX_THRESHOLD !AC70_struc(n)%ac70(t)%smcmax
                    enddo
!                  call LIS_endrun()
                 enddo
              endif          
           end do
        endif        
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     XLAI    = max(LFMASS*LAPM,AC70BIOMASSmin)

     if(sla(AC70_struc(n)%ac70(t)%vegetype).ne.0) then 
        AC70_struc(n)%ac70(t)%SumWaBal%Biomass = AC70BIOMASS(t)
     endif
  enddo

end subroutine ac70_setsoilmLAI

