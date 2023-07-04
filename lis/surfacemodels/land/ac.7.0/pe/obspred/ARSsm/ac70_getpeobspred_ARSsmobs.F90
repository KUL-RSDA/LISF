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
! !ROUTINE: Ac70_getpeobspred_ARSsmobs
!  \label{Ac70_getpeobspred_ARSsmobs}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE:
subroutine Ac70_getpeobspred_ARSsmobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_soilsMod,  only : LIS_soils
  use Ac70_lsmMod, only : Ac70_struc
  use LIS_logMod,       only : LIS_verify, LIS_logunit

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: smcField
  real, pointer          :: smc(:)
!  type(ESMF_Field)       :: smstdField
!  real, pointer          :: smstd(:)
  integer                :: t
  integer                :: i
  integer                :: status

!  write(LIS_logunit,*) '[INFO] Here 1 in Ac70_getpeobspred_ARSsmobs '

  n = 1

  call ESMF_StateGet(Obj_Func,"ARS_sm",smcField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)

!  write(LIS_logunit,*) '[INFO] Here 2 in Ac70_getpeobspred_ARSsmobs '

!  call ESMF_StateGet(Obj_Func,"ARSsm standard deviation of soil moisture",smstdField,rc=status)
!  call LIS_verify(status)
!
!  call ESMF_FieldGet(smstdField,localDE=0,farrayPtr=smstd,rc=status)
!  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc(t) = Ac70_struc(n)%ac70(t)%smc(1)
!     smstd(t) = Ac70_struc(n)%ac70(t)%smc_std
  enddo

!  write(LIS_logunit,*) '[INFO] Finished Ac70_getpeobspred_ARSsmobs, smc = ', smc

end subroutine Ac70_getpeobspred_ARSsmobs



