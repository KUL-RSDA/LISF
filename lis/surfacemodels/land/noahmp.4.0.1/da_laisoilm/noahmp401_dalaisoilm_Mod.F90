!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module noahmp401_dalaisoilm_Mod
!BOP
!
! !MODULE: noahmp401_dalaisoilm_Mod
!
! !DESCRIPTION: Updates LAI and SSM at the same time.
!  
! !REVISION HISTORY:
! 15 Dec 2018: Mahdi Navari, Sujay Kumar ; Modified for noahmp401 !
! 19 Jan 2022: Samuel Scherrer: Modified from da_soilm for da_laisoilm

! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp401_dalaisoilm_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  !public :: noahmp401_dalaism_struc
!EOP

! type, public :: dalaism_dec
!     real,    allocatable       :: model_xrange(:,:,:)
!     real,    allocatable       :: model_cdf(:,:,:)
!     real,    allocatable       :: model_mu(:)
!
!     integer                :: nbins
!     integer                :: ntimes
!     integer                :: scal
!
!  end type dalaism_dec
!  
!  type(dalaism_dec), allocatable :: noahmp401_dalaism_struc(:)

contains
!BOP
! 
! !ROUTINE: noahmp401_dalaisoilm_init
! \label{noahmp401_dalaisoilm_init}
! 
! !INTERFACE:
  subroutine noahmp401_dalaisoilm_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
!    integer                :: n 
!    character*100          :: modelcdffile(LIS_rc%nnest)
!    integer                :: status
!    integer                :: ngrid
!
!    if(.not.allocated(noahmp401_dalaism_struc)) then 
!       allocate(noahmp401_dalaism_struc(LIS_rc%nnest))
!    endif

  end subroutine noahmp401_dalaisoilm_init
end module noahmp401_dalaisoilm_Mod
