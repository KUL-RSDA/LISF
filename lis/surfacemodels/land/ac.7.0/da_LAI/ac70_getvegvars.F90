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
! !ROUTINE: ac70_getvegvars
! \label{ac70_getvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_getvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the related state prognostic variables for
!  veg data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: laiField,lfmassField

  integer                :: t
  integer                :: status
  real, pointer          :: lai(:),lfmass(:)
 
  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

!  call ESMF_StateGet(LSM_State,"LeafMass",lfmassField,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGet(lfmassField,localDE=0,farrayPtr=lfmass,rc=status)
!  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lai(t) = AC70_struc(n)%ac70(t)%lai
!     lfmass(t) = AC70_struc(n)%ac70(t)%lfmass
  enddo
  
end subroutine ac70_getvegvars

