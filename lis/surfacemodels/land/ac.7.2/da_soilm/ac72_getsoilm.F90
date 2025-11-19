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
! !ROUTINE: ac72_getsoilm
! \label{ac72_getsoilm}
!
! !REVISION HISTORY:
! 19 Nov 2025: Michel Bechtold; initial implementation
!
! !INTERFACE:
subroutine ac72_getsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture and Biomass related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  real, pointer          :: soilm1(:)
  type(ESMF_Field)       :: sm2Field
  real, pointer          :: soilm2(:)
  real, pointer       :: soilm3(:)
  real, pointer       :: soilm4(:)
  real, pointer       :: soilm5(:)
  real, pointer       :: soilm6(:)
  real, pointer       :: soilm7(:)
  real, pointer       :: soilm8(:)
  real, pointer       :: soilm9(:)
  real, pointer       :: soilm10(:)
  type(ESMF_Field)       :: sm3Field, sm4Field, sm5Field, sm6Field, sm7Field, sm8Field, sm9Field, sm10Field
  integer                :: t
  integer                :: status
  character*100          :: lsm_state_objs(2)

  ! Layer 1
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in ac72_getsoilm')
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = AC72_struc(n)%ac72(t)%smc(1)
  enddo

  ! Layer 2
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in ac72_getsoilm')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm2(t) = AC72_struc(n)%ac72(t)%smc(2)
  enddo

  ! Layer 3
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in ac72_getsoilm')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm3(t) = AC72_struc(n)%ac72(t)%smc(3)
  enddo

  ! Layer 4
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in ac72_getsoilm')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm4 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm4(t) = AC72_struc(n)%ac72(t)%smc(4)
  enddo

  ! Layer 5
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 5",sm5Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm5 in ac72_getsoilm')
  call ESMF_FieldGet(sm5Field,localDE=0,farrayPtr=soilm5,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm5 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm5(t) = AC72_struc(n)%ac72(t)%smc(5)
  enddo

  ! Layer 6
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 6",sm6Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm6 in ac72_getsoilm')
  call ESMF_FieldGet(sm6Field,localDE=0,farrayPtr=soilm6,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm6 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm6(t) = AC72_struc(n)%ac72(t)%smc(6)
  enddo

  ! Layer 7
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 7",sm7Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm7 in ac72_getsoilm')
  call ESMF_FieldGet(sm7Field,localDE=0,farrayPtr=soilm7,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm7 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm7(t) = AC72_struc(n)%ac72(t)%smc(7)
  enddo

  ! Layer 8
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 8",sm8Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm8 in ac72_getsoilm')
  call ESMF_FieldGet(sm8Field,localDE=0,farrayPtr=soilm8,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm8 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm8(t) = AC72_struc(n)%ac72(t)%smc(8)
  enddo

  ! Layer 9
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 9",sm9Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm9 in ac72_getsoilm')
  call ESMF_FieldGet(sm9Field,localDE=0,farrayPtr=soilm9,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm9 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm9(t) = AC72_struc(n)%ac72(t)%smc(9)
  enddo

  ! Layer 10
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 10",sm10Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm10 in ac72_getsoilm')
  call ESMF_FieldGet(sm10Field,localDE=0,farrayPtr=soilm10,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm10 in ac72_getsoilm')
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm10(t) = AC72_struc(n)%ac72(t)%smc(10)
  enddo

end subroutine ac72_getsoilm

