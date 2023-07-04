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
! !ROUTINE: ac70_qc_Sig0obs
! \label{ac70_qc_Sig0obs}
!
! !REVISION HISTORY:
! 26/03/2021 Sara Modanesi: Initial specifications
! 02/04/2021 Sara Modanesi: added specifications for both Sig0VV and Sig0VH S1 obs and removed flag for veg. cover
!
! !INTERFACE:
subroutine ac70_qc_Sig0VVVHobs(n,k,OBS_State)
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
  type(ESMF_Field)         :: s_vvField,s_vhField

  real, pointer            :: s_vv(:),s_vh(:)
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
  real                      :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))
  ! MB
  !real                     :: TMIN_ac(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: TMIN_ac_obs(LIS_rc%obs_ngrid(k))
  ! TMIN_ac [degC]
  type(ESMF_Field)  :: TMIN_ac_Field
  real, pointer     :: TMIN_ac(:)

!-----this part is derived from ./lis/dataassim/obs/s1_sigma/read_S1_sigma.F90
  call ESMF_StateGet(OBS_State,"Observation01",s_vvField,&
       rc=status) !
  call LIS_verify(status,&
       "ESMF_StateGet failed in ac70_qc_Sig0obs s_vv")

  call ESMF_StateGet(OBS_State,"Observation02",s_vhField,&
       rc=status) !
  call LIS_verify(status,&
       "ESMF_StateGet failed in ac70_qc_Sig0obs s_vh")

  call ESMF_FieldGet(s_vvField,localDE=0,farrayPtr=s_vv,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in ac70_qc_Sigobs s_vv")

  call ESMF_FieldGet(s_vhField,localDE=0,farrayPtr=s_vh,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in ac70_qc_Sigobs s_vh")

!---------------------------------------------------------------------------  
  ! Get TMIN_ac
  call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_TMIN_AC%varname(1)), TMIN_ac_Field, rc=status)
  call LIS_verify(status, "Ac70_f2t: error getting TMIN_ac")

  call ESMF_FieldGet(TMIN_ac_Field, localDE = 0, farrayPtr = TMIN_ac, rc = status)
  call LIS_verify(status, "Ac70_f2t: error retrieving TMIN_ac")

  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       TMIN_ac,&
       TMIN_ac_obs)


  do t = 1,LIS_rc%obs_ngrid(k)
!------------------start loop considering one obs--------------------------
     if((s_vv(t).ne.LIS_rc%udef) .and. (s_vh(t).ne.LIS_rc%udef)) then 
! MB: check for frozen soil
        if(TMIN_ac_obs(t) .lt. 1.0) then 
           s_vv(t) = LIS_rc%udef
           s_vh(t) = LIS_rc%udef
       endif
     endif
  enddo

end subroutine ac70_qc_Sig0VVVHobs

