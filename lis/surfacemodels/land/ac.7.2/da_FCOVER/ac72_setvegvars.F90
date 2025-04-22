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
! !ROUTINE: ac72_setvegvars
! \label{ac72_setvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac72_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use ac72_lsmMod
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to AquaCrop's
!  model space. The state vector consists of veg
! 
!EOP

  type(ESMF_Field)       :: AC72CCiprevField

  integer                :: t
  integer                :: status
  real, pointer          :: AC72CCiprev(:)
 
  call ESMF_StateGet(LSM_State,"AC72 CCiprev",AC72CCiprevField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(AC72CCiprevField,localDE=0,farrayPtr=AC72CCiprev,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      AC72_struc(n)%ac72(t)%CCiprev = AC72CCiprev(t)
  enddo
  
end subroutine ac72_setvegvars