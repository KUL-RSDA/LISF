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
! !ROUTINE: ac72_updatesoilm
!  \label{ac72_updatesoilm}
!
! !REVISION HISTORY:
! 19 Nov 2025: Michel Bechtold; initial implementation
!
! !INTERFACE:
subroutine ac72_updatesoilm(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to ac72
!  model space. 
! 
!EOP

  ! Declarations for soil moisture layers 1 to 10
  type(ESMF_Field)       :: sm1Field, sm1IncrField
  real, pointer          :: soilm1(:), soilmIncr1(:)

  type(ESMF_Field)       :: sm2Field, sm2IncrField
  real, pointer          :: soilm2(:), soilmIncr2(:)

  type(ESMF_Field)       :: sm3Field, sm3IncrField
  real, pointer          :: soilm3(:), soilmIncr3(:)

  type(ESMF_Field)       :: sm4Field, sm4IncrField
  real, pointer          :: soilm4(:), soilmIncr4(:)

  type(ESMF_Field)       :: sm5Field, sm5IncrField
  real, pointer          :: soilm5(:), soilmIncr5(:)

  type(ESMF_Field)       :: sm6Field, sm6IncrField
  real, pointer          :: soilm6(:), soilmIncr6(:)

  type(ESMF_Field)       :: sm7Field, sm7IncrField
  real, pointer          :: soilm7(:), soilmIncr7(:)

  type(ESMF_Field)       :: sm8Field, sm8IncrField
  real, pointer          :: soilm8(:), soilmIncr8(:)

  type(ESMF_Field)       :: sm9Field, sm9IncrField
  real, pointer          :: soilm9(:), soilmIncr9(:)

  type(ESMF_Field)       :: sm10Field, sm10IncrField
  real, pointer          :: soilm10(:), soilmIncr10(:)

  integer                :: t,i,m,gid
  integer                :: status

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  ! -------- Layer 1 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 1", sm1Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 1 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm1Field, localDE=0, farrayPtr=soilm1, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 1", sm1IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 1 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm1IncrField, localDE=0, farrayPtr=soilmIncr1, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo

  ! -------- Layer 2 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 2", sm2Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 2 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm2Field, localDE=0, farrayPtr=soilm2, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 2 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 2", sm2IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 2 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm2IncrField, localDE=0, farrayPtr=soilmIncr2, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 2 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm2(t) = soilm2(t) + soilmIncr2(t)
  enddo

  ! -------- Layer 3 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 3", sm3Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 3 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm3Field, localDE=0, farrayPtr=soilm3, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 3 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 3", sm3IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 3 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm3IncrField, localDE=0, farrayPtr=soilmIncr3, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 3 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm3(t) = soilm3(t) + soilmIncr3(t)
  enddo

  ! -------- Layer 4 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 4", sm4Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 4 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm4Field, localDE=0, farrayPtr=soilm4, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 4 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 4", sm4IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 4 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm4IncrField, localDE=0, farrayPtr=soilmIncr4, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 4 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm4(t) = soilm4(t) + soilmIncr4(t)
  enddo

  ! -------- Layer 5 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 5", sm5Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 5 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm5Field, localDE=0, farrayPtr=soilm5, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 5 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 5", sm5IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 5 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm5IncrField, localDE=0, farrayPtr=soilmIncr5, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 5 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm5(t) = soilm5(t) + soilmIncr5(t)
  enddo

  ! -------- Layer 6 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 6", sm6Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 6 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm6Field, localDE=0, farrayPtr=soilm6, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 6 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 6", sm6IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 6 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm6IncrField, localDE=0, farrayPtr=soilmIncr6, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 6 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm6(t) = soilm6(t) + soilmIncr6(t)
  enddo

  ! -------- Layer 7 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 7", sm7Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 7 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm7Field, localDE=0, farrayPtr=soilm7, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 7 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 7", sm7IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 7 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm7IncrField, localDE=0, farrayPtr=soilmIncr7, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 7 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm7(t) = soilm7(t) + soilmIncr7(t)
  enddo

  ! -------- Layer 8 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 8", sm8Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 8 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm8Field, localDE=0, farrayPtr=soilm8, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 8 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 8", sm8IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 8 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm8IncrField, localDE=0, farrayPtr=soilmIncr8, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 8 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm8(t) = soilm8(t) + soilmIncr8(t)
  enddo

  ! -------- Layer 9 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 9", sm9Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 9 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm9Field, localDE=0, farrayPtr=soilm9, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 9 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 9", sm9IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 9 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm9IncrField, localDE=0, farrayPtr=soilmIncr9, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 9 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm9(t) = soilm9(t) + soilmIncr9(t)
  enddo

  ! -------- Layer 10 --------
  call ESMF_StateGet(LSM_State, "Soil Moisture Layer 10", sm10Field, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 10 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm10Field, localDE=0, farrayPtr=soilm10, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 10 failed in ac72_updatesoilm")

  call ESMF_StateGet(LSM_Incr_State, "Soil Moisture Layer 10", sm10IncrField, rc=status)
  call LIS_verify(status, "ESMF_StateGet: Soil Moisture Layer 10 failed in ac72_updatesoilm")
  call ESMF_FieldGet(sm10IncrField, localDE=0, farrayPtr=soilmIncr10, rc=status)
  call LIS_verify(status, "ESMF_FieldGet: Soil Moisture Layer 10 failed in ac72_updatesoilm")

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     soilm10(t) = soilm10(t) + soilmIncr10(t)
  enddo

end subroutine ac72_updatesoilm

