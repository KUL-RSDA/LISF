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
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac72 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! !INTERFACE:
subroutine ac72_getsoilm(n, LSM_State)

  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use ac72_lsmMod

  implicit none

  integer, intent(in)        :: n
  type(ESMF_State)           :: LSM_State

  ! Local variables
  integer                    :: t, l, status
  integer, parameter         :: n_layers = 10  ! Replace with variable if needed elsewhere
  type(ESMF_Field)           :: smField
  real, pointer              :: soilm(:)
  character(len=100)         :: field_name

  do l = 1, n_layers
     write(field_name, '(A,I0)') "Soil Moisture Layer ", l

     call ESMF_StateGet(LSM_State, field_name, smField, rc=status)
     call LIS_verify(status, 'ESMF_StateGet failed for '//trim(field_name)//' in ac72_getsoilm')

     call ESMF_FieldGet(smField, localDE=0, farrayPtr=soilm, rc=status)
     call LIS_verify(status, 'ESMF_FieldGet failed for '//trim(field_name)//' in ac72_getsoilm')

     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        soilm(t) = AC72_struc(n)%ac72(t)%smc(l)
     end do
  end do

end subroutine ac72_getsoilm

