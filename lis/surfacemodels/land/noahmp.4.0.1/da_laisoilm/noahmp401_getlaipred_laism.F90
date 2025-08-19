!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getlaipred_laism
! \label{noahmp401_getlaipred_laism}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 19 Jan 2022: Samuel Scherrer: Modified from da_LAI for da_laisoilm
!
! !INTERFACE:
subroutine noahmp401_getlaipred_laism(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dalaisoilm_Mod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the LAI obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6)
  character*50           :: units_tp(6)
  real                   :: lai(LIS_rc%npatch(n,LIS_rc%lsm_index))


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     lai(t) = noahmp401_struc(n)%noahmp401(t)%lai
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       lai,&
       obs_pred)
  
end subroutine noahmp401_getlaipred_laism

