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
! !ROUTINE: ac72_map_veg
! \label{ac72_map_veg}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac72_map_veg(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding 
!  variables in the LSM state for veg data assimilation.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: AC72CCiprevIncrField
  type(ESMF_Field)         :: obs_FCOVER_field
  real, pointer            :: AC72CCiprevincr(:)
  real, pointer            :: FCOVERobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  character*100,allocatable    :: obs_state_objs(:)
  real, allocatable            :: ac72_CCiprev(:)

  allocate(ac72_CCiprev(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"AC72 CCiprev",AC72CCiprevIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(AC72CCiprevIncrField,localDE=0,farrayPtr=AC72CCiprevincr,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))
  
  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_FCOVER_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_FCOVER_field,localDE=0,farrayPtr=FCOVERobs,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ac72_CCiprev(t)  = AC72_struc(n)%ac72(t)%CCiprev
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(FCOVERobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index).ge.0) then 
        AC72CCiprevincr(t) = FCOVERobs(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)&
             - ac72_CCiprev(t)
     else
        AC72CCiprevincr(t) = 0 
     endif
  enddo
!  stop
  deallocate(obs_state_objs)
  deallocate(ac72_CCiprev)

end subroutine ac72_map_veg