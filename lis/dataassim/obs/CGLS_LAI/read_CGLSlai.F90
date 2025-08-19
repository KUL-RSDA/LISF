!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! !ROUTINE: read_CGLSlai
! \label{read_CGLSlai}
!
! !REVISION HISTORY:
!  03 Nov 2021    Samuel Scherrer; initial reader based on MCD152AH reader
!
! !INTERFACE: 
subroutine read_CGLSlai(n, k, OBS_State, OBS_Pert_State)
    ! !USES: 
    use ESMF
    use LIS_mpiMod
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
    use LIS_dataAssimMod
    use LIS_DAobservationsMod
    use map_utils
    use LIS_pluginIndices
    use CGLSlai_Mod, only : CGLSlai_struc

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: k
    type(ESMF_State)    :: OBS_State
    type(ESMF_State)    :: OBS_Pert_State
    !
    ! !DESCRIPTION:
    !  
    !  reads the CGLS LAI observations from NETCDF files.

    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[n] index of the nest
    !  \item[k] number of observation state
    !  \item[OBS\_State] observations state
    !  \item[OBS\_Pert\_State] observation perturbations state
    !  \end{description}
    !
    !EOP
    real, parameter        ::  minssdev = 0.05
    real, parameter        ::  maxssdev = 0.11
    real, allocatable      :: ssdev(:)
    real,  parameter       :: MAX_LAI_VALUE=10.0, MIN_LAI_VALUE=0.0001
    integer                :: status
    integer                :: grid_index
    character(len=255)     :: laiobsdir
    character(len=255)     :: fname
    integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
    real                   :: wt1, wt2,ts
    integer                :: count
    real                   :: cgmt
    real*8                 :: time
    logical                :: alarmCheck, file_exists,dataCheck
    integer                :: t,c,r,i,j,p,jj
    real,          pointer :: obsl(:)
    type(ESMF_Field)       :: laifield, pertField
    integer                :: gid(LIS_rc%obs_ngrid(k))
    integer                :: assimflag(LIS_rc%obs_ngrid(k))
    logical                :: data_update
    logical                :: data_upd_flag(LIS_npes)
    logical                :: data_upd_flag_local
    logical                :: data_upd
    real                   :: laiobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                :: fnd
    real                   :: timenow
    integer                :: prev_month
    real                   :: ndays
    real                   :: days(12)
    data days /31,28,31,30,31,30,31,31,30,31,30,31/  


    call ESMF_AttributeGet(OBS_State,"Data Directory",&
         laiobsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(OBS_State,"Data Update Status",&
         data_update, rc=status)
    call LIS_verify(status)

    data_upd = .false. 

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "CGLS LAI read alarm")

    if(alarmCheck.or.CGLSlai_struc(n)%startMode) then 
        CGLSlai_struc(n)%startMode = .false.

        call create_CGLS_LAI_filename(&
            CGLSlai_struc(n)%isresampled, CGLSlai_struc(n)%spatialres,&
            laiobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

        inquire(file=fname,exist=file_exists)          
        if(file_exists) then 
            write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
            call read_CGLS_LAI_data(n,k, fname,laiobs)
            fnd = 1
        else
            fnd = 0 
            write(LIS_logunit,*) '[WARN] Missing LAI file: ',trim(fname)
        endif
    else
        fnd = 0 
        laiobs = LIS_rc%udef
    endif

    dataCheck = .false.
    if(alarmCheck) then 
        if(fnd.eq.1) then 
            dataCheck = .true. 
        endif
    else
        fnd = 0 
        dataCheck = .false.
    endif

    if(dataCheck) then 

        call ESMF_StateGet(OBS_State,"Observation01",laifield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        call ESMF_FieldGet(laifield,localDE=0,farrayPtr=obsl,rc=status)
        call LIS_verify(status, 'Error: FieldGet')

        ! if (LIS_rc%da == 1) then
        !     ! set the data averaging factor for EnKS to the number of days
        !     ! since the last observation
        !     prev_month = mod((LIS_rc%mo - 1) + 1, 12) + 1
        !     ndays = days(prev_month) - 20.0  ! update always on the 20th
        !  else
        !     ndays = 10.0
        !  endif

        ndays = days(LIS_rc%mo)
        call ESMF_AttributeSet(OBS_State,&
            name="Data averaging factor",&
            value=ndays,rc=status)
        call LIS_verify(status)


        !-------------------------------------------------------------------------
        !  Transform data to the LSM climatology using a CDF-scaling approach
        !-------------------------------------------------------------------------     
        
        if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then

            call LIS_rescale_with_CDF_matching(     &
                 n,k,                               & 
                 CGLSlai_struc(n)%nbins,         & 
                 CGLSlai_struc(n)%ntimes,        & 
                 MAX_LAI_VALUE,                      & 
                 MIN_LAI_VALUE,                      & 
                 CGLSlai_struc(n)%model_xrange,  &
                 CGLSlai_struc(n)%obs_xrange,    &
                 CGLSlai_struc(n)%model_cdf,     &
                 CGLSlai_struc(n)%obs_cdf,       &
                 laiobs)
        endif

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         laiobs(c+(r-1)*LIS_rc%obs_lnc(k))
                endif
            enddo
        enddo

        if(fnd.eq.0) then 
            data_upd_flag_local = .false. 
        else
            data_upd_flag_local = .true. 
        endif

#if (defined SPMD)
        call MPI_ALLGATHER(data_upd_flag_local,1, &
             MPI_LOGICAL, data_upd_flag(:),&
             1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
        data_upd = .false.
        do p=1,LIS_npes
            data_upd = data_upd.or.data_upd_flag(p)
        enddo

        if(data_upd) then 

            do t=1,LIS_rc%obs_ngrid(k)
                gid(t) = t
                if(obsl(t).ne.-9999.0) then 
                    assimflag(t) = 1
                else
                    assimflag(t) = 0
                endif
            enddo

            call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .true. , rc=status)
            call LIS_verify(status)

            if(LIS_rc%obs_ngrid(k).gt.0) then 
                call ESMF_AttributeSet(laifield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(laifield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

            endif
            if(LIS_rc%dascaloption(k).eq."CDF matching") then 
                if(CGLSlai_struc(n)%useSsdevScal.eq.1) then
                    call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                         rc=status)
                    call LIS_verify(status, 'Error: StateGet Observation01')

                    allocate(ssdev(LIS_rc%obs_ngrid(k)))
                    ssdev = CGLSlai_struc(n)%ssdev_inp 

                    if(CGLSlai_struc(n)%ntimes.eq.1) then 
                        jj = 1
                    else
                        jj = LIS_rc%mo
                    endif
                    do t=1,LIS_rc%obs_ngrid(k)
                        if(CGLSlai_struc(n)%obs_sigma(t,jj).gt.0) then 
                            ssdev(t) = ssdev(t)*CGLSlai_struc(n)%model_sigma(t,jj)/&
                                 CGLSlai_struc(n)%obs_sigma(t,jj)
                            if(ssdev(t).lt.minssdev) then 
                                ssdev(t) = minssdev
                            endif
                        endif
                    enddo

                    if(LIS_rc%obs_ngrid(k).gt.0) then 
                        call ESMF_AttributeSet(pertField,"Standard Deviation",&
                             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                        call LIS_verify(status)
                    endif
                    deallocate(ssdev)
                endif
            endif

        else
            call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .false., rc=status)
            call LIS_verify(status)     
        endif
    else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)     
    endif
end subroutine read_CGLSlai

!BOP
! 
! !ROUTINE: read_CGLS_LAI_data
! \label{read_CGLS_LAI_data}
!
! !INTERFACE:
subroutine read_CGLS_LAI_data(n, k, fname, laiobs_ip)
    ! 
    ! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod
    use LIS_timeMgrMod
    use CGLSlai_Mod, only : CGLSlai_struc

    implicit none
    !
    ! !INPUT PARAMETERS: 
    ! 
    integer                       :: n 
    integer                       :: k
    character (len=*)             :: fname
    real                          :: laiobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

    ! !OUTPUT PARAMETERS:
    !
    !
    ! !DESCRIPTION: 
    !  This subroutine reads the CGLS LAI file and applies the data
    !  quality flags to filter the data. 
    !
    !  The arguments are: 
    !  \begin{description}
    !  \item[n]            index of the nest
    !  \item[k]            number of observation state
    !  \item[k]            number of observation state
    !  \item[fname]        name of the CGLS LAI file
    !  \item[laiobs\_ip]   CGLS LAI data processed to the LIS domain
    !  \end{description}
    !
    !
    !EOP

    integer                 :: lat_offset, lon_offset
    integer                 :: lai(CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr)
    integer                 :: flag(CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr)
    real                    :: lai_flagged(CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr)
    real                    :: lai_in(CGLSlai_struc(n)%nc*CGLSlai_struc(n)%nr)
    logical*1               :: lai_data_b(CGLSlai_struc(n)%nc*CGLSlai_struc(n)%nr)
    logical*1               :: laiobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                 :: c,r,t
    integer                 :: nid
    integer                 :: laiid, flagid
    integer                 :: ios

    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                                :: numLons, numLats

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    !values
    lai_data_b = .false.

    lat_offset = 1  ! no offset
    lon_offset = 1


    if (CGLSlai_struc(n)%isresampled.eq.0) then
        ! read the data from a file and optionally apply quality flags
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        ios = nf90_inq_varid(nid, 'LAI',laiid)
        call LIS_verify(ios, 'Error nf90_inq_varid: LAI')

        ios = nf90_inq_varid(nid, 'QFLAG',flagid)
        call LIS_verify(ios, 'Error nf90_inq_varid: QFLAG')

        ios = nf90_get_var(nid, laiid, lai, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr/)) 

        call LIS_verify(ios, 'Error nf90_get_var: LAI')

        ios = nf90_get_var(nid, flagid, flag, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr/))

        call LIS_verify(ios, 'Error nf90_get_var: QFLAG')

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))

        do r=1, CGLSlai_struc(n)%nr
            do c=1, CGLSlai_struc(n)%nc

                if(CGLSlai_struc(n)%qcflag.eq.1) then !apply QC flag

                    if(lai(c,r).gt.0) then
                        if (is_valid_CGLSlai_flag(flag(c,r))) then
                            lai_flagged(c,r) =&
                                 lai(c,r)*CGLSlai_struc(n)%scale
                        else
                            lai_flagged(c,r) = LIS_rc%udef
                        endif
                    else
                        lai_flagged(c,r) = LIS_rc%udef
                    endif

                else  ! no QC flag applied                

                    if(lai(c,r).gt.0) then
                        lai_flagged(c,r) =&
                             lai(c,r)*CGLSlai_struc(n)%scale
                    else
                        lai_flagged(c,r) = LIS_rc%udef
                    endif
                endif
            end do
        end do
    else
        ! if the data has been resampled, we assume that it also has been
        ! unpacked and flagged already, so we can directly read it into
        ! lai_flagged
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        ios = nf90_inq_varid(nid, 'CGLS_LAI',laiid)
        call LIS_verify(ios, 'Error nf90_inq_varid: CGLS_LAI')

        ios = nf90_get_var(nid, laiid, lai_flagged, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSlai_struc(n)%nc,CGLSlai_struc(n)%nr/)) 

        call LIS_verify(ios, 'Error nf90_get_var: CGLS_LAI')

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))

        ! the data is already read into lai_flagged, but we have to replace
        ! NaNs/invalid values with LIS_rc%udef
        do r=1, CGLSlai_struc(n)%nr
            do c=1, CGLSlai_struc(n)%nc
                if (isnan(lai_flagged(c, r))) then
                    lai_flagged(c, r) = LIS_rc%udef
                else if (.not. (0.0 < lai_flagged(c, r) .and. lai_flagged(c, r) < 20.0)) then
                    lai_flagged(c, r) = LIS_rc%udef
                endif
            end do
        end do
    endif


    ! fill lai_in and lai_data_b, which are required further on
    do r=1, CGLSlai_struc(n)%nr
        do c=1, CGLSlai_struc(n)%nc
            lai_in(c+(r-1)*CGLSlai_struc(n)%nc) = lai_flagged(c,r)
            if(lai_flagged(c,r).ne.LIS_rc%udef) then
                lai_data_b(c+(r-1)*CGLSlai_struc(n)%nc) = .true.
            else
                lai_data_b(c+(r-1)*CGLSlai_struc(n)%nc) = .false.
            endif
        enddo
    enddo

    if(LIS_rc%obs_gridDesc(k,10).lt.CGLSlai_struc(n)%dlon) then 
        write(LIS_logunit,*) '[INFO] interpolating CGLS LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Interpolate to the LIS running domain if model has finer resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
             CGLSlai_struc(n)%nc*CGLSlai_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             CGLSlai_struc(n)%rlat,CGLSlai_struc(n)%rlon,&
             CGLSlai_struc(n)%w11,CGLSlai_struc(n)%w12,&
             CGLSlai_struc(n)%w21,CGLSlai_struc(n)%w22,&
             CGLSlai_struc(n)%n11,CGLSlai_struc(n)%n12,&
             CGLSlai_struc(n)%n21,CGLSlai_struc(n)%n22,LIS_rc%udef,ios)
     else
        write(LIS_logunit,*) '[INFO] upscaling CGLS LAI',trim(fname)
        !--------------------------------------------------------------------------
        ! Upscale to the LIS running domain if model has coarser resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call upscaleByAveraging(CGLSlai_struc(n)%nc*CGLSlai_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, CGLSlai_struc(n)%n11,&
             lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
    endif

#endif

contains

    function is_valid_CGLSlai_flag(flag) result(isvalid)
        implicit none
        integer, value :: flag
        logical :: sea, filled, no_obs, lai_invalid, climato_filled, gap_filled
        ! logical :: fapar_invalid, fcover_invalid, high_lat_correction, EBF, bare
        logical :: isvalid

        sea = (iand(flag, 1) /= 0)
        filled = (iand(flag, 4) /= 0)
        no_obs = (iand(flag, 32) /= 0)
        lai_invalid = (iand(flag, 64) /= 0)
        ! fapar_invalid = (iand(flag, 128) /= 0)
        ! fcover_invalid = (iand(flag, 256) /= 0)
        ! high_lat_correction = (iand(flag, 512) /= 0)
        ! EBF = (iand(flag, 1024) /= 0)
        ! bare = (iand(flag, 2048) /= 0)
        climato_filled = (iand(flag, 4096) /= 0)
        gap_filled = (iand(flag, 8192) /= 0)

        isvalid = .not. sea &
            .and. .not. filled &
            .and. .not. no_obs &
            .and. .not. lai_invalid &
            .and. .not. climato_filled &
            .and. .not. gap_filled

    end function is_valid_CGLSlai_flag

end subroutine read_CGLS_LAI_data

!BOP
! !ROUTINE: create_CGLSlai_filename
! \label{create_CGLSlai_filename}
! 
! !INTERFACE: 
subroutine create_CGLS_LAI_filename(isresampled, res, ndir, year, month, day, filename)
    ! !USES:   

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in)              :: isresampled
    real*8, intent(in)               :: res
    character(len=*), intent(in)     :: ndir
    integer, intent(in)              :: year, month, day
    character(len=*), intent(inout)  :: filename
    ! 
    ! !DESCRIPTION: 
    !  This subroutine creates the CGLS LAI filename
    !  based on the time and date 
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[isresampled] whether the original or the resampled files are read
    !  \item[res] resolution of the files
    !  \item[ndir] name of the CGLS LAI data directory
    !  \item[year]  current year
    !  \item[month]  current month
    !  \item[day]  current day
    !  \item[filename] Generated CGLS LAI filename
    !  \end{description}
    !
    !EOP

    if (isresampled.ne.0) then
        call create_CGLSlai_filename_from_resampled(res, ndir, year, month, day, filename)
    else
         call create_CGLSlai_filename_from_original(ndir, year, month, day, filename)
    endif

contains
    subroutine create_CGLSlai_filename_from_original(ndir, year, month, day, filename)
        implicit none
        character (len=*) :: ndir
        integer, value    :: year, month, day
        character(len=*)  :: filename

        !  The naming scheme is
        !    <prefix>_<year><month><day>0000_GLOBE_<sensor>_<version>/
        !        c_gls_<prefix2>_<year><month><day>0000_GLOBE_<sensor>_<version>.nc
        !  where
        !    <prefix> is LAI or LAI_RT6
        !    <prefix2> is LAI or LAI-RT6 (corresponding always to <version>)
        !    <sensor> is VGT or PROBAV
        !    <version> is V2.0.1 or V2.0.2
        !
        !  Based on the data access portal for 1km LAI, version 2, the following
        !  combinations are available:
        !  - until 2003-06-30:
        !      <prefix> = LAI, <sensor> = VGT, <version> = V2.0.2
        !  - from 2003-07 to 2013-12-31:
        !      <prefix> = LAI, <sensor> = VGT, <version> = V2.0.1
        !  - from 2014 to 2017-05-31:
        !      <prefix> = LAI_RT6, <sensor> = PROBAV, <version> = V2.0.2
        !  - from 2017-06 to 2020-04-30:
        !      <prefix> = LAI_RT6, <sensor> = PROBAV, <version> = V2.0.1

        character (len=7) :: prefix
        character (len=7) :: prefix2
        character (len=6) :: sensor
        character (len=6) :: version
        character (len=8) :: time

        write(unit=time, fmt='((i4.4)(i2.2)(i2.2))') year, month, day

        if (year < 2003 .or. (year == 2003 .and. month <= 6)) then
            prefix = "LAI"
            prefix2 = "LAI"
            sensor = "VGT"
            version = "V2.0.2"
        else if (year < 2014) then
            prefix = "LAI"
            prefix2 = "LAI"
            sensor = "VGT"
            version = "V2.0.1"
        else if (year < 2017 .or. (year == 2017 .and. month <= 5)) then
            prefix = "LAI_RT6"
            prefix2 = "LAI-RT6"
            sensor = "PROBAV"
            version = "V2.0.2"
        else 
            prefix = "LAI_RT6"
            prefix2 = "LAI-RT6"
            sensor = "PROBAV"
            version = "V2.0.1"
        endif

        filename = trim(ndir)//'/'//&
             trim(prefix)//'_'//trim(time)//'0000_GLOBE_'//trim(sensor)//'_'//trim(version)//'/'//&
             'c_gls_'//trim(prefix2)//'_'//trim(time)//'0000_GLOBE_'//trim(sensor)//'_'//trim(version)//'.nc'

    end subroutine create_CGLSlai_filename_from_original

    subroutine create_CGLSlai_filename_from_resampled(res, ndir, year, month, day, filename)
        implicit none
        real*8, value     :: res
        character (len=*) :: ndir
        integer, value    :: year, month, day
        character(len=*)  :: filename

        !  The naming scheme is
        !    CGLS_LAI_resampled_<res>_<year>_<month>_<day>.nc
        !  where
        !    <res> is the resolution in degrees with two decimals

        character(len=5) :: resstr
        character(len=4) :: yearstr
        character(len=2) :: monthstr
        character(len=2) :: daystr

        write(unit=resstr, fmt='(f5.2)') res
        resstr=adjustl(resstr)
        write(unit=yearstr, fmt='(i4.4)') year
        write(unit=monthstr, fmt='(i2.2)') month
        write(unit=daystr, fmt='(i2.2)') day


        filename = trim(ndir)//'/'//&
             'CGLS_LAI_resampled_'//trim(resstr)//'deg_'//yearstr//'_'//monthstr//'_'//daystr//'.nc'

    end subroutine create_CGLSlai_filename_from_resampled

end subroutine create_CGLS_LAI_filename


