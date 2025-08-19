!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: Noahmp401_sfc2vod
!  \label{Noahmp401_sfc2vod}
!
! !REVISION HISTORY:
!  4 Sep 2020: Sara Modanesi; Initial Specification
! 16 Mar 2022: Samuel Scherrer; adaption for VOD observation operator
! !INTERFACE:
subroutine noahmp401_sfc2vod(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod
  use LIS_logMod,    only : LIS_verify
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use Noahmp401_lsmMod
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE, PSISAT_TABLE, BEXP_TABLE

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState

! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the noahmp401 specific surface variables
! to the VOD observation operator. 
!
!EOP
  type(ESMF_Field)    :: laiField, rzsmfield
  real, pointer       :: lai(:), rzsm(:)
  integer             :: t,status, soiltyp, nroot, iz, sm1, sm2, sm3, sm4
  real                :: sh2o, tmp_psi, dz, ztotal, bexp, eah, esat

  call ESMF_StateGet(sfcState,"Leaf Area Index",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai, rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Root Zone Soil Moisture",rzsmfield,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzsmfield,localDE=0,farrayPtr=rzsm, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      lai(t) = noahmp401_struc(n)%noahmp401(t)%lai

      sm1 = noahmp401_struc(n)%noahmp401(t)%smc(1)
      sm2 = noahmp401_struc(n)%noahmp401(t)%smc(2)
      sm3 = noahmp401_struc(n)%noahmp401(t)%smc(3)
      sm4 = noahmp401_struc(n)%noahmp401(t)%smc(4)

      ! calculate root-zone sm
      nroot = noahmp401_struc(n)%noahmp401(t)%param%nroot
      ztotal = 0.0
      rzsm(t) = 0.0
      do iz = 1, nroot
          dz = noahmp401_struc(n)%sldpth(iz)
          ztotal = ztotal + dz
          sh2o = noahmp401_struc(n)%noahmp401(t)%sh2o(iz)
          rzsm(t) = rzsm(t) + dz * sh2o
      enddo
      if (ztotal.gt.0) then
          rzsm(t) = rzsm(t) / ztotal
      endif
  enddo

end subroutine noahmp401_sfc2vod
