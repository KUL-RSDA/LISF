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
! !ROUTINE: noahmp401_qc_VODobs
! \label{noahmp401_qc_VODobs}
!
! !REVISION HISTORY:
! 26/03/2021 Sara Modanesi: Initial specifications
! 12/05/2021 Sara Modanesi: added specifications for Sig0VV S1 obs only and removed flag for veg. cover
! 18/06/2021 Michel Bechtold: assimilation flag for urban and water tiles
! 16/03/2022 Samuel Scherrer: adapted for VOD
!
! !INTERFACE:
subroutine noahmp401_qc_VODobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE
#if defined(RTMS)
  use VODOO_Mod, only : vodoo_struc
#endif


  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
#if defined(RTMS)
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the backscatter observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: vodField

  real, pointer            :: obsl(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon

! mn
  integer                  :: SOILTYP           ! soil type index [-]
  real                     :: rtm_valid(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: smc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: sh2o4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: SMCMAX(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: SMCWLT(LIS_rc%npatch(n,LIS_rc%lsm_index))

  real                     :: rtm_valid_obs(LIS_rc%obs_ngrid(k))
  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
  real                     :: sca_obs(LIS_rc%obs_ngrid(k))
!  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k)) !commented for now
  real                     :: t1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcwlt_obs(LIS_rc%obs_ngrid(k))
  real                     :: smcmax_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: smc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o1_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o2_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o3_obs(LIS_rc%obs_ngrid(k))
  real                     :: sh2o4_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))

!-----this part is derived from ./lis/dataassim/obs/s1_sigma/read_S1_sigma.F90
  call ESMF_StateGet(OBS_State,"Observation01",vodField,&
       rc=status) !
  call LIS_verify(status,&
       "ESMF_StateGet failed in noahmp401_qc_VODobs")

  call ESMF_FieldGet(vodField,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in noahmp401_qc_VODobs")

!---------------------------------------------------------------------------  

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (vodoo_struc(n)%isvalid(n,t)) then
         ! setting to a real because conversion to obs space is only defined
         ! for reals
         rtm_valid(t) = 1.0
     else
         rtm_valid(t) = 0.0
     endif
     smc1(t) = noahmp401_struc(n)%noahmp401(t)%smc(1)
     smc2(t) = noahmp401_struc(n)%noahmp401(t)%smc(2)
     smc3(t) = noahmp401_struc(n)%noahmp401(t)%smc(3)
     smc4(t) = noahmp401_struc(n)%noahmp401(t)%smc(4)

     sh2o1(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(1)
     sh2o2(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(2)
     sh2o3(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(3)
     sh2o4(t) = noahmp401_struc(n)%noahmp401(t)%sh2o(4)
!--------------------------------------------------------------------------------------------------------- 
       ! MN  NOTE:sstc contain soil and snow temperature first snow 
       !          temperature and then soil temeprature.
       !          But the number of snow layers changes from 0 to 3 
!---------------------------------------------------------------------------------------------------------
     stc1(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1)
     stc2(t) = noahmp401_struc(n)%noahmp401(t)%tslb(2)
     stc3(t) = noahmp401_struc(n)%noahmp401(t)%tslb(3)
     stc4(t) = noahmp401_struc(n)%noahmp401(t)%tslb(4)

     vegt(t) = noahmp401_struc(n)%noahmp401(t)%vegetype

     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype        
     SMCMAX(t) = SMCMAX_TABLE(SOILTYP) 
     SMCWLT(t) = SMCWLT_TABLE(SOILTYP)
  enddo

  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       rtm_valid, rtm_valid_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%prcp,&
       rainf_obs)  ! MN prcp is total precip 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%sneqv,&
       sneqv_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%snowc,&
       sca_obs)
!  call LIS_convertPatchSpaceToObsSpace(n,k,&  !out-commented for now
!       LIS_rc%lsm_index, &
!       noahmp401_struc(n)%noahmp401(:)%fveg,&
!       shdfac_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%tg,&
       t1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smcmax, &
       smcmax_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smcwlt,&
       smcwlt_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       smc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc2,&
       smc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc3,&
       smc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc4,&
       smc4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o1,&
       sh2o1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o2,&
       sh2o2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o3,&
       sh2o3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       sh2o4,&
       sh2o4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc1,&
       stc1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc2,&
       stc2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc3,&
       stc3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       stc4,&
       stc4_obs)

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)

   do t = 1,LIS_rc%obs_ngrid(k)
       !------------------start loop considering one obs--------------------------
       if(rtm_valid_obs(t).ne.1.0) then
           obsl(t) = LIS_rc%udef
       endif

       if(obsl(t).ne.LIS_rc%udef) then 
          ! MN: check for rain
          if(rainf_obs(t).gt.3E-6) then   ! Var name Noah36 --> rainf 
              obsl(t) = LIS_rc%udef
          ! MN: check for frozen soil
          elseif(abs(smc1_obs(t)- &
               sh2o1_obs(t)).gt.0.0001) then
             obsl(t) = LIS_rc%udef
          elseif(abs(smc2_obs(t)- &
               sh2o2_obs(t)).gt.0.0001) then
             obsl(t) = LIS_rc%udef
          elseif(abs(smc3_obs(t)- &
               sh2o3_obs(t)).gt.0.0001) then
             obsl(t) = LIS_rc%udef
          elseif(abs(smc4_obs(t)- &
               sh2o4_obs(t)).gt.0.0001) then
             obsl(t) = LIS_rc%udef
          elseif(stc1_obs(t).le.LIS_CONST_TKFRZ) then
             obsl(t) = LIS_rc%udef
          elseif(stc2_obs(t).le.LIS_CONST_TKFRZ) then
             obsl(t) = LIS_rc%udef
          elseif(stc3_obs(t).le.LIS_CONST_TKFRZ) then
             obsl(t) = LIS_rc%udef
          elseif(stc4_obs(t).le.LIS_CONST_TKFRZ) then
             obsl(t) = LIS_rc%udef 
          elseif(t1_obs(t).le.LIS_CONST_TKFRZ) then ! Var name Noah36 --> t1
             obsl(t) = LIS_rc%udef
  !        elseif(vegt_obs(t).le.4) then !forest types ! Var name Noah36 --> vegt
  !           obsl(t) = LIS_rc%udef
          elseif(vegt_obs(t).eq.13) then !urban ! Var name Noah36 --> vegt
             obsl(t) = LIS_rc%udef
          elseif(vegt_obs(t).eq.17) then !urban ! Var name Noah36 --> vegt
             obsl(t) = LIS_rc%udef
   ! MN: check for snow  
          elseif(sneqv_obs(t).gt.0.001) then 
             obsl(t) = LIS_rc%udef
          elseif(sca_obs(t).gt.0.0001) then  ! Var name sca 
             obsl(t) = LIS_rc%udef
   ! MN: check for green vegetation fraction NOTE: threshold incerased from 0.5 to 0.7 !commented out for now
   !      elseif(shdfac_obs(t).gt.0.7) then  ! var name Noah36 shdfac 12-month green veg. frac.  
   !          obsl(t) = LIS_rc%udef        
  !too close to the tails, could be due to scaling, so reject. !commented out for
  !now. It was written for soil moisture
  !        elseif(smcmax_obs(t)-obsl(t).lt.0.02) then 
  !           obsl(t) = LIS_rc%udef
  !        elseif(obsl(t) - smcwlt_obs(t).lt.0.02) then 
  !           obsl(t) = LIS_rc%udef
          endif
       endif
   enddo
#endif

end subroutine noahmp401_qc_VODobs
