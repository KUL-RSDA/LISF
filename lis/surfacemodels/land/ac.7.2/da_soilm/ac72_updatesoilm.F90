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
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari; Modified for ac72 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM updating with S1 backscatter w/ WCM
!
! !INTERFACE:
subroutine ac72_updatesoilm(n, LSM_State, LSM_Incr_State, n_layers)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac72_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)              :: n
  type(ESMF_State)                 :: LSM_State
  type(ESMF_State)                 :: LSM_Incr_State
  integer, parameter         :: n_layers = 10  ! Replace with variable if needed elsewhere

  integer                          :: nlayers
  character(len=64)               :: field_name
  integer                          :: l, t, i, m, gid, status
  type(ESMF_Field), allocatable    :: smField(:), smIncrField(:)
  real, pointer, dimension(:), allocatable :: soilm(:), soilmIncr(:)

  logical                         :: update_flag(LIS_rc%ngrid(n))
  real                            :: perc_violation(LIS_rc%ngrid(n))

  ! Allocate arrays
  allocate(smField(nlayers), smIncrField(nlayers))
  allocate(soilm(nlayers), soilmIncr(nlayers))

  ! Loop through each soil layer
  do l = 1, nlayers
     write(field_name, '(A,I0)') "Soil Moisture Layer ", l

     call ESMF_StateGet(LSM_State, trim(field_name), smField(l), rc=status)
     call LIS_verify(status, "ESMF_StateGet failed: "//trim(field_name)//" in ac72_updatesoilm")

     call ESMF_FieldGet(smField(l), localDE=0, farrayPtr=soilm(l), rc=status)
     call LIS_verify(status, "ESMF_FieldGet failed: "//trim(field_name)//" in ac72_updatesoilm")

     call ESMF_StateGet(LSM_Incr_State, trim(field_name), smIncrField(l), rc=status)
     call LIS_verify(status, "ESMF_StateGet (Incr) failed: "//trim(field_name)//" in ac72_updatesoilm")

     call ESMF_FieldGet(smIncrField(l), localDE=0, farrayPtr=soilmIncr(l), rc=status)
     call LIS_verify(status, "ESMF_FieldGet (Incr) failed: "//trim(field_name)//" in ac72_updatesoilm")

     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        soilm(l)(t) = soilm(l)(t) + soilmIncr(l)(t)
     end do
  end do

end subroutine ac72_updatesoilm

