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
! 19 Nov 2025: Michel Bechtold; initial implementation
!
! !INTERFACE:
subroutine ac72_getsmpred(n, k,obs_pred)
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
  real                   :: obs_tmp
  integer                :: i,t,m,gid,kk
  real                   :: inputs_tp(6), sm_out
  real                   :: w1, w2, w3
  character*50           :: units_tp(6)
  real                   :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))


  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     w1 = 1.0/6.0
     smc1(t) = w1*AC72_struc(n)%ac72(t)%smc(1) + & 
               w1*AC72_struc(n)%ac72(t)%smc(2) + &
               w1*AC72_struc(n)%ac72(t)%smc(3) + &
               w1*AC72_struc(n)%ac72(t)%smc(4) + &
               w1*AC72_struc(n)%ac72(t)%smc(5) + &
               w1*AC72_struc(n)%ac72(t)%smc(6)
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       obs_pred)

end subroutine ac72_getsmpred

