!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module ac72_prep_f

contains
  !BOP
  !
  ! !ROUTINE: AC72_ETo_calc
  ! \label{AC72_ETo_calc}
  !
  ! !REVISION HISTORY:
  !  04 NOV 2024, Louise Busschaert; initial implementation
  !

  ! !INTERFACE:
  subroutine AC72_ETo_calc(P, Tmax, Tmin, Tdew, ws, Rs, z, lat, eto)

    ! !USES:
    use LIS_constantsMod, only: LIS_CONST_PI, LIS_CONST_TKFRZ
    use LIS_coreMod, only: LIS_rc
    use LIS_tbotAdjustMod, only: LIS_tbotTimeUtil

    !
    ! !DESCRIPTION:
    !
    !  This routine computes the reference evapotranspiration (ETo)
    !  with the Penman-Monteith equation, following the guidelines
    !  of the FAO Irrigation and Drainage Paper No. 56.
    !
    !
    !EOP

    implicit none

    real, intent(in)      :: P ! kPa
    real, intent(in)      :: Tmax, Tmin, Tdew
    real, intent(in)      :: ws ! wind speed
    real, intent(in)      :: Rs ! radiation
    real, intent(in)      :: z ! elevation
    real, intent(in)      :: lat
    real, intent(inout)   :: eto ! returns eto [mm day-1]

    real                  :: Tmean
    real                  :: gamma, ea, es, slope
    real                  :: phi, dr, delta, omega, Ra, Rso
    real                  :: Rns, Rnl, Rn
    real                  :: julian_in, ratio

    ! Teamn (degC)
    Tmean = (Tmin + Tmax)/2.

    ! Psychrometric constant [kPa degC-1]
    gamma = 0.664742 * 0.001 * P

    ! Mean saturation vapour pressure (es) [kPa]
    es = (0.6108 * EXP((17.27 * Tmin) / (237.3 + Tmin)) &
         + 0.6108 * EXP((17.27 * Tmax) / (237.3 + Tmax))) / 2
    ! Actual vapor pressure (ea) [kPa]
    ea = 0.6108 * EXP((17.27 * Tdew) / (237.3 + Tdew))

    ! Slope of saturation vapour pressure curve [kPa degC-1]
    slope = (4098 * (0.6108 * EXP((17.27 * Tmean)/(Tmean + 237.3)))) &
         / (Tmean + 237.3)**2

    ! Extraterrestrial radiation (Ra) [MJ m-2 day-1]
    call LIS_tbotTimeUtil(julian_in,LIS_rc%yr) ! First get day of year
    phi = lat * LIS_CONST_PI/180 ! get latitude in rad
    dr = 1 + 0.033 * COS(2*LIS_CONST_PI/365 * julian_in) ! Inverse relative distance Earth-Sun
    delta = 0.409 * SIN(2*LIS_CONST_PI/365 * julian_in - 1.39) ! Solar declination
    omega = ACOS(-TAN(phi)*TAN(delta))
    Ra = (24 * 60/LIS_CONST_PI) * 0.082 * dr * (omega*SIN(phi)*SIN(delta) &
         + COS(phi)*COS(delta)*SIN(omega))
    if(Ra.le.0) then ! return 0 ETo for unrealisic numbers
       eto = 0
    endif

    ! Net radiation at crop surface (Rn) [MJ m-2 day-1]
    Rso = (0.75 + 2E-5 * z) * Ra
    Rns = (1 - 0.23) * Rs ! fixed value of 0.23 for albedo
    ratio = Rs/Rso
    if (ratio.gt.1) then
       ratio=1
    endif
    Rnl = 4.903E-9 * (((Tmax + LIS_CONST_TKFRZ)**4 + (Tmin + LIS_CONST_TKFRZ)**4) / 2.) &
         * (0.34 - 0.14 * SQRT(ea)) &
         * (1.35 * ratio - 0.35) ! Net longwave radiation [MJ m-2 day-1]
    Rn = Rns - Rnl

    ! Penman-Monteith --> reference evapotranspiration [mm day-1]
    eto = (0.408 * slope * Rn + (gamma * (900 / (Tmean + LIS_CONST_TKFRZ)) * ws * (es - ea))) &
         / (slope + gamma * (1 + 0.34 * ws))
    if(eto.le.0) then
       eto = 0 ! avoid negative values
    endif

  end subroutine AC72_ETo_calc

  !BOP
  !
  ! !ROUTINE: AC72_read_Trecord
  ! \label{AC72_read_Trecord}
  !
  ! !REVISION HISTORY:
  !  04 NOV 2024, Louise Busschaert; initial implementation
  !  22 APR 2025, Louise Busschaert; added plating criterion
  !

  ! !INTERFACE:
  subroutine ac72_read_Trecord(n)

    ! !USES:
    use AC72_lsmMod
    use ESMF
    use LIS_constantsMod, only: LIS_CONST_TKFRZ
    use LIS_coreMod,        only: LIS_rc, LIS_surface, LIS_resetTimeMgr
    use LIS_FORC_AttributesMod
    use LIS_logMod,         only: LIS_logunit, LIS_verify
    use LIS_metForcingMod,  only: LIS_get_met_forcing, LIS_perturb_forcing, LIS_FORC_State
    use LIS_PRIV_rcMod,     only: lisrcdec
    use LIS_timeMgrMod,     only: LIS_advance_timestep, LIS_is_last_step
    use gmaopert_Mod,       only: forcPert, forcpertdec


    !
    ! !DESCRIPTION:
    !
    !  This subroutine stores the mean temperatures for the ac72 simulation
    !  period required when AquaCrop for the determination of the GDD stages
    !
    !
    !EOP

    implicit none

    integer, intent(in)      :: n

    real, allocatable     :: daily_tmax_arr(:,:), daily_tmin_arr(:,:)
    real, allocatable     :: daily_pcp_arr(:,:)
    real, allocatable     :: subdaily_arr(:,:)
    type(lisrcdec)        :: LIS_rc_saved
    type(forcpertdec), allocatable :: forcPert_saved(:)
    integer               :: i, j, k, t, status, met_ts, m, tid
    integer               :: yr_start
    real, parameter :: lapse = -0.0065

    ! Near Surface Air Temperature [K]
    type(ESMF_Field)  :: tmpField
    real, pointer     :: tmp(:)

    ! Rainfall Rate [kg m-2 s-1]
    type(ESMF_Field)  :: pcpField
    real, pointer     :: pcp(:)

    external :: finalmetforc, initmetforc

    write(LIS_logunit,*) "[INFO] AC72: new simulation period, reading of temperature record..."

    ! Save current LIS_rc
    LIS_rc_saved = LIS_rc
    ! Save full forcing perturbation state
    allocate(forcPert_saved(LIS_rc%nforcepert))
    do k = 1, LIS_rc%nforcepert
      ! allocate deep copies
      allocate(forcPert_saved(k)%Forcepert_rseed, source=forcPert(n,k)%Forcepert_rseed)
      allocate(forcPert_saved(k)%Forcepert_ntrmdt, source=forcPert(n,k)%Forcepert_ntrmdt)
    end do

    ! Re-initialize met forcings
    do m=1,LIS_rc%nmetforc
       call finalmetforc(trim(LIS_rc%metforc(m))//char(0),m)
       call initmetforc(trim(LIS_rc%metforc(m))//char(0),m)
    enddo
    LIS_rc%rstflag(n) = 1

    met_ts = int(86400./LIS_rc%ts)

    allocate(daily_tmax_arr(LIS_rc%npatch(n,LIS_rc%lsm_index),366))
    allocate(daily_tmin_arr(LIS_rc%npatch(n,LIS_rc%lsm_index),366))
    allocate(subdaily_arr(LIS_rc%npatch(n,LIS_rc%lsm_index),met_ts))
    daily_tmax_arr = 0
    daily_tmin_arr = 0

    if (AC72_struc(n)%Rainfall_crit) then
        allocate(daily_pcp_arr(LIS_rc%npatch(n,LIS_rc%lsm_index),366))
        daily_pcp_arr = 0
    endif

    ! Set LIS_rc time to beginning of simulation period (in case of restart)
    ! Check in which year the simulation did start (assuming a 365-366 sim period)
    if (AC72_struc(n)%Sim_AnnualStartMonth.gt.LIS_rc%mo) then
       yr_start = LIS_rc%yr - 1
    else
       yr_start = LIS_rc%yr
    endif
    ! Set LIS_rc to reset clock
    LIS_rc%syr = yr_start
    LIS_rc%smo = AC72_struc(n)%Sim_AnnualStartMonth
    LIS_rc%sda = AC72_struc(n)%Sim_AnnualStartDay
    LIS_rc%shr = LIS_rc%hr+1
    call LIS_resetTimeMgr
    day_loop: do i=1,366
       do j=1,met_ts
          ! read met forcing
          call LIS_get_met_forcing(n)
          call LIS_perturb_forcing(n)

          ! Get Tair
          call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Tair%varname(1)), tmpField, rc=status)
          call LIS_verify(status, "AC72_prep_f: error getting Tair")

          call ESMF_FieldGet(tmpField, localDE = 0, farrayPtr = tmp, rc = status)
          call LIS_verify(status, "AC72_prep_f: error retrieving Tair")

          ! Store temperatures
          subdaily_arr(:,j) = tmp

          if (AC72_struc(n)%Rainfall_crit) then
            call LIS_perturb_forcing(n)
            ! Get and store rainfall (for sowing/planting based on rainfall criterion)
            call ESMF_StateGet(LIS_FORC_State(n), trim(LIS_FORC_Rainf%varname(1)), pcpField, rc=status)
            call LIS_verify(status, "AC72_f2t: error getting Rainf")

            call ESMF_FieldGet(pcpField, localDE = 0, farrayPtr = pcp, rc = status)
            call LIS_verify(status, "AC72_f2t: error retrieving Rainf")

            daily_pcp_arr(:,i) = daily_pcp_arr(:,i) + pcp*LIS_rc%ts
          endif

          ! Change LIS time to the next meteo time step
          call LIS_advance_timestep(LIS_rc)
       enddo
       ! Store daily max and min temperatures
       daily_tmax_arr(:,i) = maxval(subdaily_arr,2)
       daily_tmin_arr(:,i) = minval(subdaily_arr,2)
       if ((LIS_rc%da.eq.AC72_struc(n)%Sim_AnnualStartDay)&
            .and.(LIS_rc%mo.eq.AC72_struc(n)%Sim_AnnualStartMonth)&
            .and.(LIS_rc%hr.eq.LIS_rc_saved%hr+1)&
            .and.(i.ne.1)) exit day_loop
       ! Exit if we reach end of sim period
       ! but still include the last hour
    enddo day_loop

    deallocate(subdaily_arr)

    ! lapse-rate correct the temperature records
    daily_tmax_arr(:,:) = daily_tmax_arr(:,:) + (lapse * (2 -AC72_struc(n)%forchgt_tq))
    daily_tmin_arr(:,:) = daily_tmin_arr(:,:) + (lapse * (2 -AC72_struc(n)%forchgt_tq))
    ! Assign Tmax and Tmin arrays to AC72_struc
    do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
       tid = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%tile_id
       ! AquaCrop rounds meteo input to 4 decimals
       AC72_struc(n)%ac72(t)%Tmax_record = anint((daily_tmax_arr(tid,:)-LIS_CONST_TKFRZ)*10000)/10000
       AC72_struc(n)%ac72(t)%Tmin_record = anint((daily_tmin_arr(tid,:)-LIS_CONST_TKFRZ)*10000)/10000
       if (AC72_struc(n)%Rainfall_crit) then
           AC72_struc(n)%ac72(t)%pcp_record =  anint(daily_pcp_arr(tid,:)*10000)/10000
       endif
    enddo

    deallocate(daily_tmax_arr)
    deallocate(daily_tmin_arr)
    if (AC72_struc(n)%Rainfall_crit) then
        deallocate(daily_pcp_arr)
    endif

    ! Reset LIS_rc
    LIS_rc = LIS_rc_saved
    ! Restore full forcing perturbation state
    do k = 1, LIS_rc%nforcepert
      forcPert(n,k)%Forcepert_rseed  = forcPert_saved(k)%Forcepert_rseed
      forcPert(n,k)%Forcepert_ntrmdt = forcPert_saved(k)%Forcepert_ntrmdt
    end do
    deallocate(forcPert_saved)

    ! Set LIS_rc to reset clock
    LIS_rc%syr = LIS_rc_saved%yr
    LIS_rc%smo = LIS_rc_saved%mo
    LIS_rc%sda = LIS_rc_saved%da
    LIS_rc%shr = LIS_rc_saved%hr
    call LIS_resetTimeMgr
    ! Re-initialize met forcings
    do m=1,LIS_rc%nmetforc
       call finalmetforc(trim(LIS_rc%metforc(m))//char(0),m)
       call initmetforc(trim(LIS_rc%metforc(m))//char(0),m)
    enddo
    LIS_rc%rstflag(n) = 1 ! For met forcings

    ! Check if end of LIS run
    if (LIS_is_last_step(LIS_rc)) then
       LIS_rc%endtime = 1
    endif

    write(LIS_logunit,*) "[INFO] AC72: new simulation period, reading of temperature record... Done!"
  end subroutine ac72_read_Trecord

!BOP
  !
  ! !ROUTINE: AC72_search_start_Temp
  ! \label{AC72_search_start_Temp}
  !
  ! !REVISION HISTORY:
  !  22 APR 2025, Louise Busschaert; initial implementation (planting criterion)
  !
integer function ac72_search_start_Temp(startsim_julian_days, startcrop_julian_days, crit_window, &
                                        Temp_crit_tmin, Temp_crit_days, Temp_crit_occurrence, Tmin_record)

    !
    ! !DESCRIPTION:
    !
    !  This function finds the planting/sowing date based on a temperature
    !  criterion defined in the lis.config
    !
    !
    !EOP

    implicit none
    ! Input parameters
    integer, intent(in) :: startsim_julian_days, startcrop_julian_days, crit_window, &
                           Temp_crit_days, Temp_crit_occurrence
    integer, intent(in) :: Temp_crit_tmin
    real, intent(in) :: Tmin_record(366)  ! Yearly min temperature data (1-366)

    ! Local variables
    integer :: search_start, search_end, t, i, occurrence_count
    logical :: consecutive_met

    ! Define search start and end indices
    search_start = startcrop_julian_days - startsim_julian_days + 1
    search_end = search_start + crit_window - 1

    ! Initialize occurrence counter
    occurrence_count = 0
    t = search_start

    ! Loop through the search window with non-overlapping windows
    do while (t <= search_end - Temp_crit_days + 1)
        consecutive_met = .true.

        ! Check Tmin_record for Temp_crit_days consecutive days
        do i = 0, Temp_crit_days - 1
            if (Tmin_record(t + i) < Temp_crit_tmin) then
                consecutive_met = .false.
                exit
            end if
        end do

        if (consecutive_met) then
            occurrence_count = occurrence_count + 1
            if (occurrence_count == Temp_crit_occurrence) then
                ac72_search_start_Temp = startcrop_julian_days + (t + i - search_start)
                return
            end if
            ! Move `t` to the end of the window
            t = t + i
        end if

        ! Move to the next day
        t = t + 1
    end do

    ! If condition is not met, return the last possible day
    ac72_search_start_Temp = startcrop_julian_days + (search_end - search_start)

end function ac72_search_start_Temp

!BOP
  !
  ! !ROUTINE: AC72_search_start_Rainfall
  ! \label{AC72_search_start_Rainfall}
  !
  ! !REVISION HISTORY:
  !  22 APR 2025, Louise Busschaert; initial implementation (planting criterion)
  !
integer function ac72_search_start_Rainfall(startsim_julian_days, startcrop_julian_days, crit_window, &
                                            Rainfall_crit_amount, Rainfall_crit_days, Rainfall_crit_occurrence, pcp_record)

    !
    ! !DESCRIPTION:
    !
    !  This function finds the planting/sowing date based on a rainfall
    !  criterion defined in the lis.config
    !
    !
    !EOP

    implicit none
    ! Input parameters
    integer, intent(in) :: startsim_julian_days, startcrop_julian_days, crit_window, &
                           Rainfall_crit_days, Rainfall_crit_occurrence
    integer, intent(in) :: Rainfall_crit_amount
    real, intent(in) :: pcp_record(366)  ! Yearly precipitation data (1-366)

    ! Local variables
    integer :: search_start, search_end, t, i, occurrence_count
    real    :: sum_pcp

    ! Define search start and end indices
    search_start = startcrop_julian_days - startsim_julian_days + 1
    search_end = search_start + crit_window - 1

    ! Initialize occurrence counter
    occurrence_count = 0
    t = search_start

    ! Loop through the search window
    do while (t <= search_end - Rainfall_crit_days + 1)
        sum_pcp = 0.0  ! Reset sum for this window

        ! Compute sum of rainfall over Rainfall_crit_days
        do i = 0, Rainfall_crit_days - 1
            sum_pcp = sum_pcp + pcp_record(t + i)
            if (sum_pcp > Rainfall_crit_amount) then
                occurrence_count = occurrence_count + 1
                if (occurrence_count == Rainfall_crit_occurrence) then
                    ac72_search_start_Rainfall = startcrop_julian_days + (t + i - search_start)
                    return
                end if
                ! Move `t` to the end of the window
                t = t + i
                exit  ! Exit inner loop early
            end if
        end do

        ! Move to the next day
        t = t + 1
    end do

    ! If condition is not met, return the last possible day
    ac72_search_start_Rainfall = startcrop_julian_days + (search_end - search_start)

end function ac72_search_start_Rainfall

end module ac72_prep_f
