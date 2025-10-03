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
! !ROUTINE: ac72_getsmpred
! \label{ac72_getsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari 
! 23 Nov 2022: Michel Bechtold, Modified for ac72 
!
! !INTERFACE:
subroutine ac72_getsmpred(n, k, obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use ac72_lsmMod
  use ac72_dasoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  !-------------------------------
  ! Parameters
  !-------------------------------
  integer, parameter     :: n_layers = 10   ! number of soil layers to use

  !-------------------------------
  ! Local variables
  !-------------------------------
  integer                :: t, l
  real                   :: smc(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: w(n_layers)

  ! Define weights equally distributed over n_layers
  w(:) = 1.0 / real(n_layers)

  ! Compute weighted soil moisture for each patch
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     smc(t) = 0.0
     do l=1, n_layers
        smc(t) = smc(t) + w(l) * AC72_struc(n)%ac72(t)%smc(l)
     enddo
  enddo

  ! Convert to observation ensemble space
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc,&
       obs_pred)

end subroutine ac72_getsmpred

