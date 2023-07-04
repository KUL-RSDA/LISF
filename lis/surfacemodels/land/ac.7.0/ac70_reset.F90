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
!
! !ROUTINE: Ac70_reset
! \label{Ac70_reset}
!
! !REVISION HISTORY:
!  22 Feb 2018: Soni Yatheendradas; Initial version
! 
! !INTERFACE:
subroutine Ac70_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use Ac70_lsmMod
  use LIS_logMod,       only : LIS_verify, LIS_logunit

!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Ac70 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Ac70.
!  
!EOP
  implicit none
  integer                 :: tt,n
  integer                 :: status


  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)                        &
          'Ac70 resetting'

     ! initialize forcing variables to zeros
     do tt=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
         AC70_struc(n)%ac70(tt)%lwdown = 0.0
         AC70_struc(n)%ac70(tt)%swdown = 0.0
         AC70_struc(n)%ac70(tt)%psurf = 0.0
         AC70_struc(n)%ac70(tt)%prcp = 0.0
         AC70_struc(n)%ac70(tt)%tair = 0.0
         AC70_struc(n)%ac70(tt)%qair = 0.0
         AC70_struc(n)%ac70(tt)%wind_e = 0.0
         AC70_struc(n)%ac70(tt)%wind_n = 0.0
         AC70_struc(n)%ac70(tt)%PREC_ac = 0.0
         AC70_struc(n)%ac70(tt)%TMIN_ac = 0.0
         AC70_struc(n)%ac70(tt)%TMAX_ac = 0.0
         AC70_struc(n)%ac70(tt)%ETo_ac = 0.0
     enddo ! end of tile (tt) loop
     AC70_struc(n)%forc_count = 0
     
  enddo ! do n=1,LIS_rc%nnest
end subroutine Ac70_reset
