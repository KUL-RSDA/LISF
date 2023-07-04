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
! !ROUTINE: ac70_getsoilmLAI
! \label{ac70_getsoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac70 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
! !INTERFACE:
subroutine ac70_getsoilmLAI(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture and LAI related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: AC70BIOMASSField,lfmassField
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: AC70BIOMASS(:)
  character*100          :: lsm_state_objs(2)


  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in ac70_getsoilmLAI')
  call ESMF_StateGet(LSM_State,"AC70 BIOMASS",AC70BIOMASSField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for AC70BIOMASS in ac70_getsoilmLAI')

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in ac70_getsoilmLAI')
  call ESMF_FieldGet(AC70BIOMASSField,localDE=0,farrayPtr=AC70BIOMASS,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for AC70BIOMASS in ac70_getsoilmLAI')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = AC70_struc(n)%ac70(t)%smc(1)
     !CCIprev(t) = AC70_struc(n)%ac70(t)%CCiprev
     AC70BIOMASS(t) = AC70_struc(n)%ac70(t)%SumWaBal%Biomass
  enddo

end subroutine ac70_getsoilmLAI

