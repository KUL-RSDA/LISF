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
! !ROUTINE: ac70_setvegvars
! \label{ac70_setvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use ac70_lsmMod
  use AC_VEG_PARAMETERS_70
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to noah's
!  model space. The state vector consists of veg
! 
!EOP

  type(ESMF_Field)       :: AC70CCiprevField

  integer                :: t
  integer                :: status
  real, pointer          :: AC70CCiprev(:)
 
  call ESMF_StateGet(LSM_State,"AC70 CCiprev",AC70CCiprevField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(AC70CCiprevField,localDE=0,farrayPtr=AC70CCiprev,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      AC70_struc(n)%ac70(t)%CCiprev = AC70CCiprev(t)
  enddo
  
end subroutine ac70_setvegvars


