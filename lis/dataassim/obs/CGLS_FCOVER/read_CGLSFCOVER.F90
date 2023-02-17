!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"

!BOP
! !ROUTINE: read_CGLSFCOVER
! \label{read_CGLSFCOVER}
!
! !REVISION HISTORY:
!  15 Feb 2023    Zdenko Heyvaert; initial reader based on Samuel Scherrer's CGLS LAI reader
!
! !INTERFACE: 
subroutine read_CGLSFCOVER(n, k, OBS_State, OBS_Pert_State)
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
    use CGLSFCOVER_Mod, only : CGLSFCOVER_struc

    implicit none
    ! !ARGUMENTS: 
    integer, intent(in) :: n 
    integer, intent(in) :: k
    type(ESMF_State)    :: OBS_State
    type(ESMF_State)    :: OBS_Pert_State
    !
    ! !DESCRIPTION:
    !  
    !  reads the CGLS FCOVER observations from NETCDF files.

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
    real,  parameter       :: MAX_FCOVER_VALUE=1.0, MIN_FCOVER_VALUE=0.0001
    integer                :: status
    integer                :: grid_index
    character(len=255)     :: FCOVERobsdir
    character(len=255)     :: fname
    integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
    real                   :: wt1, wt2,ts
    integer                :: count
    real                   :: cgmt
    real*8                 :: time
    logical                :: alarmCheck, file_exists,dataCheck
    integer                :: t,c,r,i,j,p,jj
    real,          pointer :: obsl(:)
    type(ESMF_Field)       :: FCOVERfield, pertField
    integer                :: gid(LIS_rc%obs_ngrid(k))
    integer                :: assimflag(LIS_rc%obs_ngrid(k))
    logical                :: data_update
    logical                :: data_upd_flag(LIS_npes)
    logical                :: data_upd_flag_local
    logical                :: data_upd
    real                   :: FCOVERobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                :: fnd
    real                   :: timenow
    integer                :: prev_month
    real                   :: ndays
    real                   :: days(12)
    data days /31,28,31,30,31,30,31,31,30,31,30,31/  


    call ESMF_AttributeGet(OBS_State,"Data Directory",&
         FCOVERobsdir, rc=status)
    call LIS_verify(status)
    call ESMF_AttributeGet(OBS_State,"Data Update Status",&
         data_update, rc=status)
    call LIS_verify(status)

    data_upd = .false. 

    alarmCheck = LIS_isAlarmRinging(LIS_rc, "CGLS FCOVER read alarm")

    if(alarmCheck.or.CGLSFCOVER_struc(n)%startMode) then 
        CGLSFCOVER_struc(n)%startMode = .false.

        call create_CGLS_FCOVER_filename(&
            CGLSFCOVER_struc(n)%isresampled, CGLSFCOVER_struc(n)%spatialres,&
            FCOVERobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

        inquire(file=fname,exist=file_exists)          
        if(file_exists) then 
            write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
            call read_CGLS_FCOVER_data(n,k, fname,FCOVERobs)
            fnd = 1
        else
            fnd = 0 
            write(LIS_logunit,*) '[WARN] Missing FCOVER file: ',trim(fname)
        endif
    else
        fnd = 0 
        FCOVERobs = LIS_rc%udef
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

        call ESMF_StateGet(OBS_State,"Observation01",FCOVERfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        call ESMF_FieldGet(FCOVERfield,localDE=0,farrayPtr=obsl,rc=status)
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
                 CGLSFCOVER_struc(n)%nbins,         & 
                 CGLSFCOVER_struc(n)%ntimes,        & 
                 MAX_FCOVER_VALUE,                      & 
                 MIN_FCOVER_VALUE,                      & 
                 CGLSFCOVER_struc(n)%model_xrange,  &
                 CGLSFCOVER_struc(n)%obs_xrange,    &
                 CGLSFCOVER_struc(n)%model_cdf,     &
                 CGLSFCOVER_struc(n)%obs_cdf,       &
                 FCOVERobs)
        endif

        obsl = LIS_rc%udef 
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                    obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                         FCOVERobs(c+(r-1)*LIS_rc%obs_lnc(k))
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
                call ESMF_AttributeSet(FCOVERfield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(FCOVERfield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

            endif
            if(LIS_rc%dascaloption(k).eq."CDF matching") then 
                if(CGLSFCOVER_struc(n)%useSsdevScal.eq.1) then
                    call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                         rc=status)
                    call LIS_verify(status, 'Error: StateGet Observation01')

                    allocate(ssdev(LIS_rc%obs_ngrid(k)))
                    ssdev = CGLSFCOVER_struc(n)%ssdev_inp 

                    if(CGLSFCOVER_struc(n)%ntimes.eq.1) then 
                        jj = 1
                    else
                        jj = LIS_rc%mo
                    endif
                    do t=1,LIS_rc%obs_ngrid(k)
                        if(CGLSFCOVER_struc(n)%obs_sigma(t,jj).gt.0) then 
                            ssdev(t) = ssdev(t)*CGLSFCOVER_struc(n)%model_sigma(t,jj)/&
                                 CGLSFCOVER_struc(n)%obs_sigma(t,jj)
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
end subroutine read_CGLSFCOVER

!BOP
! 
! !ROUTINE: read_CGLS_FCOVER_data
! \label{read_CGLS_FCOVER_data}
!
! !INTERFACE:
subroutine read_CGLS_FCOVER_data(n, k, fname, FCOVERobs_ip)
    ! 
    ! !USES:   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod
    use LIS_timeMgrMod
    use CGLSFCOVER_Mod, only : CGLSFCOVER_struc

    implicit none
    !
    ! !INPUT PARAMETERS: 
    ! 
    integer                       :: n 
    integer                       :: k
    character (len=*)             :: fname
    real                          :: FCOVERobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

    ! !OUTPUT PARAMETERS:
    !
    !
    ! !DESCRIPTION: 
    !  This subroutine reads the CGLS FCOVER file and applies the data
    !  quality flags to filter the data. 
    !
    !  The arguments are: 
    !  \begin{description}
    !  \item[n]            index of the nest
    !  \item[k]            number of observation state
    !  \item[k]            number of observation state
    !  \item[fname]        name of the CGLS FCOVER file
    !  \item[FCOVERobs\_ip]   CGLS FCOVER data processed to the LIS domain
    !  \end{description}
    !
    !
    !EOP

    integer                 :: lat_offset, lon_offset
    integer                 :: FCOVER(CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr)
    integer                 :: flag(CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr)
    real                    :: FCOVER_flagged(CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr)
    real                    :: FCOVER_in(CGLSFCOVER_struc(n)%nc*CGLSFCOVER_struc(n)%nr)
    logical*1               :: FCOVER_data_b(CGLSFCOVER_struc(n)%nc*CGLSFCOVER_struc(n)%nr)
    logical*1               :: FCOVERobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
    integer                 :: c,r,t
    integer                 :: nid
    integer                 :: FCOVERid, flagid
    integer                 :: ios

    integer, dimension(nf90_max_var_dims) :: dimIDs
    integer                                :: numLons, numLats

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    !values
    FCOVER_data_b = .false.

    lat_offset = 1  ! no offset
    lon_offset = 1


    if (CGLSFCOVER_struc(n)%isresampled.eq.0) then
        ! read the data from a file and optionally apply quality flags
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        ios = nf90_inq_varid(nid, 'FCOVER',FCOVERid)
        call LIS_verify(ios, 'Error nf90_inq_varid: FCOVER')

        ios = nf90_inq_varid(nid, 'QFLAG',flagid)
        call LIS_verify(ios, 'Error nf90_inq_varid: QFLAG')

        ios = nf90_get_var(nid, FCOVERid, FCOVER, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr/)) 

        call LIS_verify(ios, 'Error nf90_get_var: FCOVER')

        ios = nf90_get_var(nid, flagid, flag, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr/))

        call LIS_verify(ios, 'Error nf90_get_var: QFLAG')

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))

        do r=1, CGLSFCOVER_struc(n)%nr
            do c=1, CGLSFCOVER_struc(n)%nc

                if(CGLSFCOVER_struc(n)%qcflag.eq.1) then !apply QC flag

                    if(FCOVER(c,r).gt.0) then
                        if (is_valid_CGLSFCOVER_flag(flag(c,r))) then
                            FCOVER_flagged(c,r) =&
                                 FCOVER(c,r)*CGLSFCOVER_struc(n)%scale
                        else
                            FCOVER_flagged(c,r) = LIS_rc%udef
                        endif
                    else
                        FCOVER_flagged(c,r) = LIS_rc%udef
                    endif

                else  ! no QC flag applied                

                    if(FCOVER(c,r).gt.0) then
                        FCOVER_flagged(c,r) =&
                             FCOVER(c,r)*CGLSFCOVER_struc(n)%scale
                    else
                        FCOVER_flagged(c,r) = LIS_rc%udef
                    endif
                endif
            end do
        end do
    else
        ! if the data has been resampled, we assume that it also has been
        ! unpacked and flagged already, so we can directly read it into
        ! FCOVER_flagged
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        ios = nf90_inq_varid(nid, 'CGLS_FCOVER',FCOVERid)
        call LIS_verify(ios, 'Error nf90_inq_varid: CGLS_FCOVER')

        ios = nf90_get_var(nid, FCOVERid, FCOVER_flagged, &
             start=(/lon_offset,lat_offset/), &
             count=(/CGLSFCOVER_struc(n)%nc,CGLSFCOVER_struc(n)%nr/)) 

        call LIS_verify(ios, 'Error nf90_get_var: CGLS_FCOVER')

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))

        ! the data is already read into FCOVER_flagged, but we have to replace
        ! NaNs/invalid values with LIS_rc%udef
        do r=1, CGLSFCOVER_struc(n)%nr
            do c=1, CGLSFCOVER_struc(n)%nc
                if (isnan(FCOVER_flagged(c, r))) then
                    FCOVER_flagged(c, r) = LIS_rc%udef
                else if (.not. (0.0 < FCOVER_flagged(c, r) .and. FCOVER_flagged(c, r) < 20.0)) then
                    FCOVER_flagged(c, r) = LIS_rc%udef
                endif
            end do
        end do
    endif


    ! fill FCOVER_in and FCOVER_data_b, which are required further on
    do r=1, CGLSFCOVER_struc(n)%nr
        do c=1, CGLSFCOVER_struc(n)%nc
            FCOVER_in(c+(r-1)*CGLSFCOVER_struc(n)%nc) = FCOVER_flagged(c,r)
            if(FCOVER_flagged(c,r).ne.LIS_rc%udef) then
                FCOVER_data_b(c+(r-1)*CGLSFCOVER_struc(n)%nc) = .true.
            else
                FCOVER_data_b(c+(r-1)*CGLSFCOVER_struc(n)%nc) = .false.
            endif
        enddo
    enddo

    if(LIS_rc%obs_gridDesc(k,10).lt.CGLSFCOVER_struc(n)%dlon) then 
        write(LIS_logunit,*) '[INFO] interpolating CGLS FCOVER',trim(fname)
        !--------------------------------------------------------------------------
        ! Interpolate to the LIS running domain if model has finer resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             FCOVER_data_b, FCOVER_in, FCOVERobs_b_ip, FCOVERobs_ip, &
             CGLSFCOVER_struc(n)%nc*CGLSFCOVER_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             CGLSFCOVER_struc(n)%rlat,CGLSFCOVER_struc(n)%rlon,&
             CGLSFCOVER_struc(n)%w11,CGLSFCOVER_struc(n)%w12,&
             CGLSFCOVER_struc(n)%w21,CGLSFCOVER_struc(n)%w22,&
             CGLSFCOVER_struc(n)%n11,CGLSFCOVER_struc(n)%n12,&
             CGLSFCOVER_struc(n)%n21,CGLSFCOVER_struc(n)%n22,LIS_rc%udef,ios)
     else
        write(LIS_logunit,*) '[INFO] upscaling CGLS FCOVER',trim(fname)
        !--------------------------------------------------------------------------
        ! Upscale to the LIS running domain if model has coarser resolution
        ! than observations
        !-------------------------------------------------------------------------- 
        call upscaleByAveraging(CGLSFCOVER_struc(n)%nc*CGLSFCOVER_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             LIS_rc%udef, CGLSFCOVER_struc(n)%n11,&
             FCOVER_data_b,FCOVER_in, FCOVERobs_b_ip, FCOVERobs_ip)
    endif

#endif

contains

    function is_valid_CGLSFCOVER_flag(flag) result(isvalid)
        implicit none
        integer, value :: flag
        logical :: sea, filled, no_obs, FCOVER_invalid, climato_filled, gap_filled
        ! logical :: fapar_invalid, FCOVER_invalid, high_lat_correction, EBF, bare
        logical :: isvalid

        sea = (iand(flag, 1) /= 0)
        filled = (iand(flag, 4) /= 0)
        no_obs = (iand(flag, 32) /= 0)
        FCOVER_invalid = (iand(flag, 64) /= 0)
        ! fapar_invalid = (iand(flag, 128) /= 0)
        ! FCOVER_invalid = (iand(flag, 256) /= 0)
        ! high_lat_correction = (iand(flag, 512) /= 0)
        ! EBF = (iand(flag, 1024) /= 0)
        ! bare = (iand(flag, 2048) /= 0)
        climato_filled = (iand(flag, 4096) /= 0)
        gap_filled = (iand(flag, 8192) /= 0)

        isvalid = .not. sea &
            .and. .not. filled &
            .and. .not. no_obs &
            .and. .not. FCOVER_invalid &
            .and. .not. climato_filled &
            .and. .not. gap_filled

    end function is_valid_CGLSFCOVER_flag

end subroutine read_CGLS_FCOVER_data

!BOP
! !ROUTINE: create_CGLSFCOVER_filename
! \label{create_CGLSFCOVER_filename}
! 
! !INTERFACE: 
subroutine create_CGLS_FCOVER_filename(isresampled, res, ndir, year, month, day, filename)
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
    !  This subroutine creates the CGLS FCOVER filename
    !  based on the time and date 
    ! 
    !  The arguments are: 
    !  \begin{description}
    !  \item[isresampled] whether the original or the resampled files are read
    !  \item[res] resolution of the files
    !  \item[ndir] name of the CGLS FCOVER data directory
    !  \item[year]  current year
    !  \item[month]  current month
    !  \item[day]  current day
    !  \item[filename] Generated CGLS FCOVER filename
    !  \end{description}
    !
    !EOP

    if (isresampled.ne.0) then
        call create_CGLSFCOVER_filename_from_resampled(res, ndir, year, month, day, filename)
    else
         call create_CGLSFCOVER_filename_from_original(ndir, year, month, day, filename)
    endif

contains
    subroutine create_CGLSFCOVER_filename_from_original(ndir, year, month, day, filename)
        implicit none
        character (len=*) :: ndir
        integer, value    :: year, month, day
        character(len=*)  :: filename

        !  The naming scheme is
        !    <prefix>_<year><month><day>0000_GLOBE_<sensor>_<version>/
        !        c_gls_<prefix2>_<year><month><day>0000_GLOBE_<sensor>_<version>.nc
        !  where
        !    <prefix> is FCOVER or FCOVER_RT6
        !    <prefix2> is FCOVER or FCOVER-RT6 (corresponding always to <version>)
        !    <sensor> is VGT or PROBAV
        !    <version> is V2.0.1 or V2.0.2
        !
        !  Based on the data access portal for 1km FCOVER, version 2, the following
        !  combinations are available:
        !  - until 2003-06-30:
        !      <prefix> = FCOVER, <sensor> = VGT, <version> = V2.0.2
        !  - from 2003-07 to 2013-12-31:
        !      <prefix> = FCOVER, <sensor> = VGT, <version> = V2.0.1
        !  - from 2014 to 2017-05-31:
        !      <prefix> = FCOVER_RT6, <sensor> = PROBAV, <version> = V2.0.2
        !  - from 2017-06 to 2020-04-30:
        !      <prefix> = FCOVER_RT6, <sensor> = PROBAV, <version> = V2.0.1

        character (len=7) :: prefix
        character (len=7) :: prefix2
        character (len=6) :: sensor
        character (len=6) :: version
        character (len=8) :: time

        write(unit=time, fmt='((i4.4)(i2.2)(i2.2))') year, month, day

        if (year < 2003 .or. (year == 2003 .and. month <= 6)) then
            prefix = "FCOVER"
            prefix2 = "FCOVER"
            sensor = "VGT"
            version = "V2.0.2"
        else if (year < 2014) then
            prefix = "FCOVER"
            prefix2 = "FCOVER"
            sensor = "VGT"
            version = "V2.0.1"
        else if (year < 2017 .or. (year == 2017 .and. month <= 5)) then
            prefix = "FCOVER_RT6"
            prefix2 = "FCOVER-RT6"
            sensor = "PROBAV"
            version = "V2.0.2"
        else 
            prefix = "FCOVER_RT6"
            prefix2 = "FCOVER-RT6"
            sensor = "PROBAV"
            version = "V2.0.1"
        endif

        filename = trim(ndir)//'/'//&
             trim(prefix)//'_'//trim(time)//'0000_GLOBE_'//trim(sensor)//'_'//trim(version)//'/'//&
             'c_gls_'//trim(prefix2)//'_'//trim(time)//'0000_GLOBE_'//trim(sensor)//'_'//trim(version)//'.nc'

    end subroutine create_CGLSFCOVER_filename_from_original

    subroutine create_CGLSFCOVER_filename_from_resampled(res, ndir, year, month, day, filename)
        implicit none
        real*8, value     :: res
        character (len=*) :: ndir
        integer, value    :: year, month, day
        character(len=*)  :: filename

        !  The naming scheme is
        !    CGLS_FCOVER_resampled_<res>_<year>_<month>_<day>.nc
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
             'CGLS_FCOVER_resampled_'//trim(resstr)//'.deg_'//yearstr//'_'//monthstr//'_'//daystr//'.nc'

    end subroutine create_CGLSFCOVER_filename_from_resampled

end subroutine create_CGLS_FCOVER_filename


