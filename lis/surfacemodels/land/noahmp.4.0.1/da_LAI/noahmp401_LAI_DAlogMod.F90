!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! 16Feb12  Ben Zaitchik; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1
! 07 Feb 2022: Samuel Scherrer; adapted from tws_DAlogMod

module noahmp401_LAI_DAlogMod
  
  use LIS_constantsMod,  only : LIS_CONST_RHOFW
  use ESMF
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------
  public :: noahmp401_LAI_DAlog
!-----------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------
  public :: NOAHMPpred_struc
!EOP

  type, public ::NOAHMPpred_dec
     
     real,allocatable ::lai(:,:)
     
  end type NOAHMPpred_dec
  
  type (NOAHMPpred_dec),allocatable :: NOAHMPpred_struc(:)
  
contains 
  
  subroutine noahmp401_LAI_DAlog(n)
    
    ! USES:
    use LIS_coreMod, only : LIS_rc,LIS_surface
    use LIS_timeMgrMod
    use noahmp401_lsmMod
    use LIS_logMod, only : LIS_logunit, LIS_verify
    !      use smootherDA_runMod, only : smootherDA_increments_mode
    implicit none
      
    ! ARGUMENTS:  
    integer, intent(in)    :: n 
      
    ! DESCRIPTION:
    ! Calculates LAI at the times when CGLS LAI observations are available
    
    integer                  :: i,m,t,gid,d          
    integer                  :: yr,mo,da,hr,mn,ss
    integer                  :: yr1, mo1, da1
    integer                  :: yr2, mo2, da2
    integer                  :: yr3, mo3, da3
    integer                  :: tw_tmp1, tw_tmp2
    type(ESMF_Time)          :: tTime1,tTime2,tTime3
    type(ESMF_TimeInterval)  :: tw1, tw2
    integer                  :: status
    real                     :: days(12)
    data days /31,28,31,30,31,30,31,31,30,31,30,31/  

    if(LIS_rc%DAincrMode(n).eq.0) then
       if(LIS_rc%twInterval.eq.2592000.0) then 
          ! if the time interval was set to 1 month, we log LAI at the 10th, the
          ! 20th and the last day of a month, as in CGLS LAI
          if((LIS_rc%da.eq.1).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then
             call noahmp401_LAI_reset_log_vars(n)
          end if
          
          mo = LIS_rc%mo
          if(((LIS_rc%da.eq.10).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%da.eq.20).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%da.eq.days(mo)).and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0))) then
             
             d = nint((LIS_rc%da)/10.0)
             call noahmp401_LAI_set_log_vars(n, d)
          endif
       
       else
          ! time interval set to something else than 1 month, 3 logging times
          ! are derived from 
          call ESMF_TimeGet(LIS_twMidTime, yy = yr, &
               mm = mo, &
               dd = da, &
               h  = hr, &
               m  = mn,& 
               s  = ss, & 
               calendar = LIS_calendar, & 
               rc = status)

          if((LIS_rc%mo.eq.mo).and.(LIS_rc%da.eq.da) &
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then
              call noahmp401_LAI_reset_log_vars(n)
          end if
          
          tw_tmp1 = nint(LIS_rc%twInterval/3.0)
          tw_tmp2 = tw_tmp1/2
          call ESMF_TimeIntervalSet(tw1,s=tw_tmp1,rc=status)
          call ESMF_TimeIntervalSet(tw2,s=tw_tmp2,rc=status)
          
          tTime1 = LIS_twMidTime + tw2
          tTime2 = tTime1 + tw1
          tTime3 = tTime2 + tw1
          
          call ESMF_TimeGet(tTime1,yy=yr1,mm=mo1,dd=da1,calendar=LIS_calendar,&
               rc=status)
          call ESMF_TimeGet(tTime2,yy=yr2,mm=mo2,dd=da2,calendar=LIS_calendar,&
               rc=status)
          call ESMF_TimeGet(tTime3,yy=yr3,mm=mo3,dd=da3,calendar=LIS_calendar,&
               rc=status)
          
          if(&
               ((LIS_rc%yr.eq.yr1).and.(LIS_rc%mo.eq.mo1).and.(LIS_rc%da.eq.da1)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%yr.eq.yr2).and.(LIS_rc%mo.eq.mo2).and.(LIS_rc%da.eq.da2)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)).or. &
               ((LIS_rc%yr.eq.yr3).and.(LIS_rc%mo.eq.mo3).and.(LIS_rc%da.eq.da3)&
               .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) & 
               ) then 
               
             d = -1
             if((LIS_rc%yr.eq.yr1).and.(LIS_rc%mo.eq.mo1).and.(LIS_rc%da.eq.da1)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 1
             elseif((LIS_rc%yr.eq.yr2).and.(LIS_rc%mo.eq.mo2).and.(LIS_rc%da.eq.da2)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 2
             elseif((LIS_rc%yr.eq.yr3).and.(LIS_rc%mo.eq.mo3).and.(LIS_rc%da.eq.da3)&
                  .and.(LIS_rc%hr.eq.12).and.(LIS_rc%mn.eq.0)) then 
                d = 3
             endif

             call noahmp401_LAI_set_log_vars(n, d)
          endif
          
       endif
    endif

  end subroutine noahmp401_LAI_DAlog

  subroutine noahmp401_LAI_reset_log_vars(n)
    ! USES:
    use LIS_coreMod, only : LIS_rc,LIS_surface
    use LIS_timeMgrMod
    use noahmp401_lsmMod
    use LIS_logMod, only : LIS_logunit, LIS_verify

    implicit none

    ! ARGUMENTS
    integer, intent(in)    :: n

    ! DESCRIPTION:
    ! Allocates and resets the variable array used for logging in
    ! NOAHMPpred_struc

    if(.not.allocated(NOAHMPpred_struc)) then 
        allocate(NOAHMPpred_struc(LIS_rc%nnest))
        allocate(NOAHMPpred_struc(n)%lai(3,&
             LIS_rc%npatch(n,LIS_rc%lsm_index)))
    endif
    NOAHMPpred_struc(n)%lai = 0.0

  end subroutine noahmp401_LAI_reset_log_vars

  subroutine noahmp401_LAI_set_log_vars(n, d)
    ! USES:
    use LIS_coreMod, only : LIS_rc,LIS_surface
    use LIS_timeMgrMod
    use noahmp401_lsmMod
    use LIS_logMod, only : LIS_logunit, LIS_verify

    implicit none

    ! ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(in)    :: d

    ! DESCRIPTION:
    ! Sets the variable to be logged in NOAHMPpred_struc

    integer :: t

    write(LIS_logunit,*)'[INFO] logging obspred data for CGLS LAI DA'

    NOAHMPpred_struc(n)%lai(d,:) = 0.0
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       NOAHMPpred_struc(n)%lai(d,t) = NOAHMP401_struc(n)%noahmp401(t)%lai
    enddo

  end subroutine noahmp401_LAI_set_log_vars
  
end module noahmp401_LAI_DAlogMod
