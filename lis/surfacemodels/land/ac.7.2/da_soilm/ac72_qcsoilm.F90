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
! !ROUTINE: ac72_qcsoilm
! \label{ac72_qcsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac72 
!
! !INTERFACE:
subroutine ac72_qcsoilm(n, LSM_State)

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac72_lsmMod

  implicit none

  integer, intent(in)        :: n
  type(ESMF_State)           :: LSM_State

  integer                    :: t, l, status
  integer, parameter         :: n_layers = 10  ! Make this dynamic if needed
  type(ESMF_Field)           :: smField
  real, pointer              :: soilm(:)
  real                       :: smmax(n_layers), smmin(n_layers)
  character(len=100)         :: field_name

  integer                    :: gid
  logical                    :: update_flag(LIS_rc%ngrid(n))
  real                       :: perc_violation(LIS_rc%ngrid(n))
  integer                    :: N_ens
  real                       :: state_tmp(LIS_rc%nensem(n)), state_mean

  !--- Soil Moisture QC ---
  do l = 1, n_layers
     write(field_name, '(A,I0)') "Soil Moisture Layer ", l

     call ESMF_StateGet(LSM_State, field_name, smField, rc=status)
     call LIS_verify(status, "ESMF_StateGet for "//trim(field_name)//" failed in ac72_qcsoilm")

     call ESMF_FieldGet(smField, localDE=0, farrayPtr=soilm, rc=status)
     call LIS_verify(status, "ESMF_FieldGet for "//trim(field_name)//" failed in ac72_qcsoilm")

     call ESMF_AttributeGet(smField, "Max Value", smmax(l), rc=status)
     call LIS_verify(status, "ESMF_AttributeGet Max Value failed for "//trim(field_name))

     call ESMF_AttributeGet(smField, "Min Value", smmin(l), rc=status)
     call LIS_verify(status, "ESMF_AttributeGet Min Value failed for "//trim(field_name))

     do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        if (soilm(t) > smmax(l)) soilm(t) = smmax(l)
        if (soilm(t) < smmin(l)) soilm(t) = smmin(l)
     end do
  end do

end subroutine ac72_qcsoilm
