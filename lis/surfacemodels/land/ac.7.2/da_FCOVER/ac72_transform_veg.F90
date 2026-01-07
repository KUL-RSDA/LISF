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
! !ROUTINE: ac72_transform_veg
! \label{ac72_transform_veg}
!
! !REVISION HISTORY:
!  22 Dec 2017: Sujay Kumar, Initial specification 
!
! !INTERFACE:
subroutine ac72_transform_veg(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use ac72_lsmMod
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the veg state
!  to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

!nothing to be done

#if 0 
  type(ESMF_Field)         :: obs_FCOVER_field
  real, pointer            :: FCOVERobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status
  
  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  call LIS_verify(status, 'attributeget error in ac72_transform_veg')
  call ESMF_StateGet(OBS_State,"Observation01",obs_FCOVER_field,&
       rc=status)
  call LIS_verify(status,'stateget error in ac72_transform_veg')
  call ESMF_FieldGet(obs_FCOVER_field,localDE=0,farrayPtr=FCOVERobs,rc=status)
  call LIS_verify(status,'fieldget error in ac72_transform_veg')

  ! If using 8th mesh LAI data, convert it from inches to meters.
  ! 16th mesh LAI data are already in meters.
  do t=1,N_obs_size
     if(FCOVERobs(t).ne.LIS_rc%udef) then 
        FCOVERobs(t) = FCOVERobs(t)
     endif
  enddo
#endif
end subroutine ac72_transform_veg