!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_merra2_elev
! \label{read_merra2_elev}
!
! !REVISION HISTORY:
!
!  22 Sep 2016: James Geiger; Initial Specification
!
! !INTERFACE:
subroutine read_merra2_elev(n,findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_fileIOMod,     only : LIS_read_param
  use LIS_logMod,        only : LIS_logunit

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

! !DESCRIPTION:
!
!  Opens and reads MERRA2 model elevation to the LIS grid.
!  The data will be used to perform any topographical adjustments
!  to the forcing.
!
!  The elevation file needs to be preprocessed to fit the running 
!  resolution and domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[LIS\_readData](\ref{LIS_readData}) \newline
!    Abstract method to read the elevation of the forcing
!    data in the same map projection used in LIS.
!  \end{description}
!EOP

  integer :: c,r
  real    :: go(LIS_rc%lnc(n),LIS_rc%lnr(n))

!  if ( trim(LIS_rc%met_ecor(findex)) .ne. "none") then 
  if ( LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
        LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" .or. &
        LIS_rc%met_ecor(findex) == "micromet" ) then

     write(LIS_logunit,*) '[INFO] Reading the MERRA2 elevation map '
     
     call LIS_read_param(n,"ELEV_MERRA2",go)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_merra2_elev
