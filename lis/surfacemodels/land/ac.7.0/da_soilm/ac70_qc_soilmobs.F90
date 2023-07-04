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
! !ROUTINE: ac70_qc_soilmobs
! \label{ac70_qc_soilmobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 1 Aug 2016: Mahdi Navari; Modified for Ac70 
!
! !INTERFACE:
subroutine ac70_qc_soilmobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod   
  use ac70_lsmMod
  use module_sf_noahaclsm_36  !, only: MAXSMC !MN


  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_sm_field

  real, pointer            :: smobs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon

! mn
  integer                 :: SOILTYP           ! soil type index [-]
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

  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
  real                     :: sca_obs(LIS_rc%obs_ngrid(k))
  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k))
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
  ! MB
  !real                     :: TMIN_ac(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: TMIN_ac_obs(LIS_rc%obs_ngrid(k))
  ! TMIN_ac [degC]
  type(ESMF_Field)  :: TMIN_ac_Field
  real, pointer     :: TMIN_ac(:)

  call ESMF_StateGet(OBS_State,"Observation01",obs_sm_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in ac70_qc_soilmobs")
  call ESMF_FieldGet(obs_sm_field,localDE=0,farrayPtr=smobs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in ac70_qc_soilmobs")
  
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc1(t) = ac70_struc(n)%ac70(t)%smc(1)

     SOILTYP = AC70_struc(n)%ac70(t)%soiltype        
     SMCMAX(t)  = MAXSMC (SOILTYP) 
     SMCWLT(t) = WLTSMC (SOILTYP)
  enddo

  ! Get TMIN_ac
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_TMIN_AC%varname(1)), TMIN_ac_Field, rc=status)
  call LIS_verify(status, "Ac70_f2t: error getting TMIN_ac")

  call ESMF_FieldGet(TMIN_ac_Field, localDE = 0, farrayPtr = TMIN_ac, rc = status)
  call LIS_verify(status, "Ac70_f2t: error retrieving TMIN_ac")

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
       TMIN_ac,&
       TMIN_ac_obs)

  do t = 1,LIS_rc%obs_ngrid(k)

     if(smobs(t).ne.LIS_rc%udef) then 
! MB: check for frozen soil
        if(TMIN_ac_obs(t) .lt. 1.0) then 
           smobs(t) = LIS_rc%udef
!too close to the tails, could be due to scaling, so reject. 
        elseif(smcmax_obs(t)-smobs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
        elseif(smobs(t) - smcwlt_obs(t).lt.0.02) then 
           smobs(t) = LIS_rc%udef
        endif
     endif
  enddo

end subroutine ac70_qc_soilmobs

