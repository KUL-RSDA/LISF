!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SMAP_AC72rzmc
! \label{read_SMAP_AC72rzmc}
!
! !REVISION HISTORY:
!  17 Nov 2025: Michel Bechtold, Initial Specification
!
! !INTERFACE:
subroutine read_SMAP_AC72rzmc(n,k,  OBS_State, OBS_Pert_State)
! !USES:
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use SMAP_AC72_rzmc_Mod, only : SMAP_AC72_rzmc_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  Reader for SMAP RZMC product optimized for AquaCrop 
!EOP
  real, parameter        ::  minssdev = 0.01
  real, parameter        ::  maxssdev = 0.11
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir, fname
  logical :: alarmCheck, file_exists
  integer :: t,c,r,p,jj
  real, pointer :: obsl(:)
  type(ESMF_Field) :: smfield, pertfield
  integer :: gid(LIS_rc%obs_ngrid(k)), assimflag(LIS_rc%obs_ngrid(k))
  real    :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical :: data_update, data_upd_flag(LIS_npes), data_upd_flag_local, data_upd
  real    :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real    :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real    :: dt, lon, lhour, gmt
  integer :: zone, fnd

  call ESMF_AttributeGet(OBS_State,"Data Directory", smobsdir, rc=status)
  call LIS_verify(status, 'ESMF_AttributeGet failed in read_SMAP_AC72rzmc')
  call ESMF_AttributeGet(OBS_State,"Data Update Status", data_update, rc=status)
  call LIS_verify(status, 'ESMF_AttributeGet failed in read_SMAP_AC72rzmc')

  data_upd = .false.; obs_unsc = LIS_rc%udef
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMAP_AC72 read alarm")

  if (alarmCheck .or. SMAP_AC72_rzmc_struc(n)%startMode) then
     SMAP_AC72_rzmc_struc(n)%startMode = .false.
     SMAP_AC72_rzmc_struc(n)%smobs = LIS_rc%udef
     smobs  = LIS_rc%udef
     SMAP_AC72_rzmc_struc(n)%smtime = -1

     call create_SMAP_AC72rzmc_filename(smobsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)
     inquire(file=fname, exist=file_exists)
     if (file_exists) then
        write(LIS_logunit,*) 'Reading ', trim(fname)
        call read_SMAP_AC72_data(n,k,fname,smobs)
     else
        write(LIS_logunit,*) 'WARNING: No SMAP RZMC file found for date: ', LIS_rc%yr, LIS_rc%mo, LIS_rc%da
     endif

     SMAP_AC72_rzmc_struc(n)%smobs  = LIS_rc%udef
     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1) then
              if (smobs(c+(r-1)*LIS_rc%obs_lnc(k)) .gt. 0.0) then
                 SMAP_AC72_rzmc_struc(n)%smobs(c,r) = smobs(c+(r-1)*LIS_rc%obs_lnc(k))
                 lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                 lhour = 12.0
                 call LIS_localtime2gmt (gmt,lon,lhour,zone)
                 SMAP_AC72_rzmc_struc(n)%smtime(c,r) = gmt
              endif
           endif
        enddo
     enddo
  endif

  call ESMF_StateGet(OBS_State,"Observation01",smfield, rc=status); call LIS_verify(status, 'StateGet Observation01')
  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status);   call LIS_verify(status, 'FieldGet')

  fnd = 0; sm_current = LIS_rc%udef
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1) then
           if (SMAP_AC72_rzmc_struc(n)%midnight_assimilation) then
              dt = (LIS_rc%gmt)*3600.0
           else
              dt = (LIS_rc%gmt - SMAP_AC72_rzmc_struc(n)%smtime(c,r))*3600.0
           endif
           if ((dt .ge. 0.0 .and. dt .lt. LIS_rc%ts) .or. &
               (LIS_rc%lsm .eq. "AquaCrop.7.2")) then
              sm_current(c,r) = SMAP_AC72_rzmc_struc(n)%smobs(c,r)
              obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = sm_current(c,r)
              if (sm_current(c,r) .ne. LIS_rc%udef) fnd = 1
           endif
        endif
     enddo
  enddo

  if (SMAP_AC72_rzmc_struc(n)%ntimes.gt.1 .and. SMAP_AC72_rzmc_struc(n)%cdf_read_opt.eq.1) then
   if (.not. SMAP_AC72_rzmc_struc(n)%cdf_read_mon .or. &
       LIS_rc%da.eq.1 .and. LIS_rc%hr.eq.0 .and. LIS_rc%mn.eq.0 .and. LIS_rc%ss.eq.0) then

      call LIS_readMeanSigmaData(n,k, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
            SMAP_AC72_rzmc_struc(n)%modelcdffile, "SoilMoist", &
            SMAP_AC72_rzmc_struc(n)%model_mu, SMAP_AC72_rzmc_struc(n)%model_sigma, LIS_rc%mo)

      call LIS_readMeanSigmaData(n,k, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
            SMAP_AC72_rzmc_struc(n)%obscdffile, "SoilMoist", &
            SMAP_AC72_rzmc_struc(n)%obs_mu, SMAP_AC72_rzmc_struc(n)%obs_sigma, LIS_rc%mo)

      call LIS_readCDFdata(n,k, SMAP_AC72_rzmc_struc(n)%nbins, SMAP_AC72_rzmc_struc(n)%ntimes, &
            LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%modelcdffile, "SoilMoist", &
            SMAP_AC72_rzmc_struc(n)%model_xrange, SMAP_AC72_rzmc_struc(n)%model_cdf, LIS_rc%mo)

      call LIS_readCDFdata(n,k, SMAP_AC72_rzmc_struc(n)%nbins, SMAP_AC72_rzmc_struc(n)%ntimes, &
            LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%obscdffile, "SoilMoist", &
            SMAP_AC72_rzmc_struc(n)%obs_xrange, SMAP_AC72_rzmc_struc(n)%obs_cdf, LIS_rc%mo)

      SMAP_AC72_rzmc_struc(n)%cdf_read_mon = .true.
   endif
  endif

  if (trim(LIS_rc%dascaloption(k)).eq."CDF matching" .and. fnd.ne.0) then
   if (SMAP_AC72_rzmc_struc(n)%ntimes.gt.1 .and. SMAP_AC72_rzmc_struc(n)%cdf_read_opt.eq.1) then
      call LIS_rescale_with_CDF_matching( n,k, SMAP_AC72_rzmc_struc(n)%nbins, 1, &
           MAX_SM_VALUE, MIN_SM_VALUE, &
           SMAP_AC72_rzmc_struc(n)%model_xrange, SMAP_AC72_rzmc_struc(n)%obs_xrange, &
           SMAP_AC72_rzmc_struc(n)%model_cdf, SMAP_AC72_rzmc_struc(n)%obs_cdf, sm_current)
   else
      call LIS_rescale_with_CDF_matching( n,k, SMAP_AC72_rzmc_struc(n)%nbins, &
           SMAP_AC72_rzmc_struc(n)%ntimes, MAX_SM_VALUE, MIN_SM_VALUE, &
           SMAP_AC72_rzmc_struc(n)%model_xrange, SMAP_AC72_rzmc_struc(n)%obs_xrange, &
           SMAP_AC72_rzmc_struc(n)%model_cdf, SMAP_AC72_rzmc_struc(n)%obs_cdf, sm_current)
   endif
  endif

  obsl = LIS_rc%udef
  do r=1, LIS_rc%obs_lnr(k)
     do c=1, LIS_rc%obs_lnc(k)
        if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1) then
           obsl(LIS_obs_domain(n,k)%gindex(c,r)) = sm_current(c,r)
        endif
     enddo
  enddo

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"//trim(LIS_SMAP_AC72rzmcobsId)//char(0), n, k, OBS_state)
  call LIS_checkForValidObs(n,k,obsl,fnd,sm_current)

  if (fnd.eq.0) then
     data_upd_flag_local = .false.
  else
     data_upd_flag_local = .true.
  endif

#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local,1, MPI_LOGICAL, data_upd_flag(:), 1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
  data_upd = .false.; do p=1,LIS_npes; data_upd = data_upd .or. data_upd_flag(p); enddo

  if (data_upd .and. alarmCheck) then
     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if (obsl(t) .ne. -9999.0) then; assimflag(t) = 1; else; assimflag(t) = 0; endif
     enddo
     call ESMF_AttributeSet(OBS_State,"Data Update Status", .true., rc=status); call LIS_verify(status)
     if (LIS_rc%obs_ngrid(k).gt.0) then
        call ESMF_AttributeSet(smField,"Grid Number", gid, itemCount=LIS_rc%obs_ngrid(k), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(smField,"Assimilation Flag", assimflag, itemCount=LIS_rc%obs_ngrid(k), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(smfield,"Unscaled Obs", obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status); call LIS_verify(status)
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status", .false., rc=status); call LIS_verify(status)
  endif
end subroutine read_SMAP_AC72rzmc

!BOP
! !ROUTINE: read_SMAP_AC72_data
! \label{read_SMAP_AC72_data}
subroutine read_SMAP_AC72_data(n, k, fname, smobs_ip)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod,          only : LIS_rc
  use LIS_logMod,           only : LIS_verify, LIS_logunit
  use SMAP_AC72_rzmc_Mod,   only : SMAP_AC72_rzmc_struc
  implicit none
  integer, intent(in)           :: n, k
  character (len=*), intent(in) :: fname
  real, intent(out)             :: smobs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  ! NetCDF ids / status
  integer :: nid, sm_varid, lat_dimid, lon_dimid, time_dimid, lat_varid, lon_varid
  integer :: ios

  ! Variable shape / positions
  integer :: vrank, vdimids(4), dimlen, d, dimid
  integer :: nlat, nlon, nt
  integer :: pos_lat, pos_lon, pos_time
  integer :: dlen3(3), t_index
  integer :: start3(3), count3(3)

  ! Data buffers
  real, allocatable :: latv(:), lonv(:)
  real, allocatable :: sm2(:,:)          ! 2-D slice (lat,lon)
  real, allocatable :: var3(:,:,:)       ! full 3-D in file order
  real, allocatable :: src_data(:)
  logical*1, allocatable :: src_mask(:)
  logical*1, allocatable :: dst_mask(:)

  ! Attributes / QC
  real :: fillv, scale, offset
  logical :: has_fill_nan
  logical :: has_scale, has_offset, is_nan

  ! Grid/build-map vars
  integer :: r, c
  real :: dlat, dlon, lat0, lon0, latN, lonN
  real :: gridDesci(50)
  logical :: need_build_map

  ! Logging helpers
  real, allocatable :: row(:)
  integer :: n_valid_src, n_valid_dst
  real :: smin, smax, dmin, dmax

  ! Target grid longitude wrap detection (optional)
  real :: tgt_lon0, tgt_lonN
  logical :: need_lon_wrap

  !----------------------
  ! Open and basic inquiry
  !----------------------
  ios = nf90_open(trim(fname), NF90_NOWRITE, nid)
  call LIS_verify(ios,'Error opening file '//trim(fname))

  ! Canonical dimension and variable ids (lat/lon/time may be missing in legacy)
  ios = nf90_inq_dimid(nid, 'lat',  lat_dimid);   if (ios /= nf90_noerr) lat_dimid  = -1
  ios = nf90_inq_dimid(nid, 'lon',  lon_dimid);   if (ios /= nf90_noerr) lon_dimid  = -1
  ios = nf90_inq_dimid(nid, 'time', time_dimid);  if (ios /= nf90_noerr) time_dimid = -1

  ios = nf90_inq_varid(nid, 'sm',  sm_varid);  call LIS_verify(ios, 'Error nf90_inq_varid: sm')
  ios = nf90_inq_varid(nid, 'lat', lat_varid); call LIS_verify(ios, 'Error nf90_inq_varid: lat')
  ios = nf90_inq_varid(nid, 'lon', lon_varid); call LIS_verify(ios, 'Error nf90_inq_varid: lon')

  ! Rank & dim ids attached to 'sm'
  ios = nf90_inquire_variable(nid, sm_varid, ndims=vrank, dimids=vdimids)
  call LIS_verify(ios, 'Error nf90_inquire_variable: sm')

  ! Defaults
  nlat = -1; nlon = -1; nt = 1
  pos_lat  = 0; pos_lon  = 0; pos_time = 0

  ! Sizes and positions by id match
  do d = 1, vrank
     dimid = vdimids(d)
     ios = nf90_inquire_dimension(nid, dimid, len=dimlen)
     call LIS_verify(ios, 'Error inquiring sm dim length')
     if (dimid == lat_dimid) then
        nlat = dimlen; pos_lat = d
     else if (dimid == lon_dimid) then
        nlon = dimlen; pos_lon = d
     else if (dimid == time_dimid) then
        nt = dimlen;   pos_time = d
     end if
  end do

  ! Fallback: assume last two dims are spatial; infer time as the remaining one (3D only)
  if ((nlat < 0 .or. nlon < 0) .and. vrank >= 2) then
     ios = nf90_inquire_dimension(nid, vdimids(vrank-1), len=nlat); call LIS_verify(ios, 'Error inferring nlat')
     ios = nf90_inquire_dimension(nid, vdimids(vrank  ), len=nlon); call LIS_verify(ios, 'Error inferring nlon')
     pos_lat = vrank-1
     pos_lon = vrank
     if (vrank == 3 .and. pos_time == 0) then
        pos_time = 6 - pos_lat - pos_lon   ! positions {1,2,3}
        ios = nf90_inquire_dimension(nid, vdimids(pos_time), len=nt); call LIS_verify(ios, 'Error inferring nt')
     end if
  end if

  ! Coordinates
  allocate(latv(nlat), lonv(nlon))
  ios = nf90_get_var(nid, lat_varid, latv); call LIS_verify(ios,'Error nf90_get_var: lat')
  ios = nf90_get_var(nid, lon_varid, lonv); call LIS_verify(ios,'Error nf90_get_var: lon')

  ! Attributes
  fillv     = -9999.0
  scale     = 1.0
  offset    = 0.0
  has_scale  = (nf90_get_att(nid, sm_varid, 'scale_factor', scale) == nf90_noerr)
  has_offset = (nf90_get_att(nid, sm_varid, 'add_offset',  offset) == nf90_noerr)
  ios = nf90_get_att(nid, sm_varid, '_FillValue', fillv)  
  has_fill_nan = (ios == nf90_noerr) .and. (fillv .ne. fillv)   ! true if NaN

  !----------------------
  ! Read sm into a 2-D (lat,lon) slice
  !----------------------
  allocate(sm2(nlat,nlon))

  if (nt < 1) then
     smobs_ip = LIS_rc%udef
     ios = nf90_close(nid); call LIS_verify(ios,'Error closing file '//trim(fname))
     return
  end if

  ! Query lengths in the variable’s order and read full 3-D
  do d = 1, 3
     ios = nf90_inquire_dimension(nid, vdimids(d), len=dlen3(d))
     call LIS_verify(ios, 'Error inquiring sm dim length (dlen3)')
  end do
  allocate(var3(dlen3(1), dlen3(2), dlen3(3)))
  ios = nf90_get_var(nid, sm_varid, var3)
  call LIS_verify(ios, 'Error nf90_get_var: sm (3D full read)')

  t_index = 1

  ! Copy first time slice into sm2(lat,lon)
  if (pos_time == 1) then                ! (time, lat, lon)
     do r=1,nlat; do c=1,nlon; sm2(r,c) = var3(t_index, r, c); end do; end do
  else if (pos_time == 2) then           ! (lat, time, lon)
     do r=1,nlat; do c=1,nlon; sm2(r,c) = var3(r, t_index, c); end do; end do
  else if (pos_time == 3) then           ! (*, *, time) — need pos_lat/pos_lon
     if (pos_lat==1 .and. pos_lon==2) then
        do r=1,nlat; do c=1,nlon; sm2(r,c) = var3(r, c, t_index); end do; end do
     else if (pos_lat==2 .and. pos_lon==1) then
        do r=1,nlat; do c=1,nlon; sm2(r,c) = var3(c, r, t_index); end do; end do
     else
        call LIS_verify(-1, 'Unexpected (lat,lon) positions with time last')
     end if
  else
     call LIS_verify(-1, 'Could not determine time dimension position')
  end if

  deallocate(var3)


  ! Done with file
  ios = nf90_close(nid); call LIS_verify(ios,'Error closing file '//trim(fname))

  !----------------------
  ! Optional scaling and masking prep
  !----------------------
  if (has_scale .or. has_offset) then
     do r=1,nlat
        do c=1,nlon
           if (sm2(r,c) /= fillv) sm2(r,c) = sm2(r,c)*scale + offset
        end do
     end do
  end if

  ! Optional: wrap source longitudes to 0..360 if target grid is 0..360
  tgt_lon0 = LIS_rc%obs_gridDesc(k,5)   ! convention: lon start
  tgt_lonN = LIS_rc%obs_gridDesc(k,8)   ! lon end
  need_lon_wrap = (tgt_lon0 >= 0.0 .and. tgt_lonN > 180.0)
  if (need_lon_wrap) then
     do c=1,nlon
        if (lonv(c) < 0.0) lonv(c) = lonv(c) + 360.0
     end do
  end if

  !----------------------
  ! Build/refresh source grid mapping if needed
  !----------------------
  need_build_map = .false.
  if (.not. allocated(SMAP_AC72_rzmc_struc(n)%rlat)) need_build_map = .true.
  if (allocated(SMAP_AC72_rzmc_struc(n)%rlat)) then
     if (size(SMAP_AC72_rzmc_struc(n)%rlat) /= nlat*nlon) need_build_map = .true.
  end if

  if (need_build_map) then
     SMAP_AC72_rzmc_struc(n)%ecvnr = nlat
     SMAP_AC72_rzmc_struc(n)%ecvnc = nlon

     if (nlat > 1) then; dlat = abs(latv(2) - latv(1)); else; dlat = 0.0; endif
     if (nlon > 1) then; dlon = abs(lonv(2) - lonv(1)); else; dlon = 0.0; endif
     lat0 = latv(1); lon0 = lonv(1)
     latN = latv(nlat); lonN = lonv(nlon)

     if (allocated(SMAP_AC72_rzmc_struc(n)%rlat)) deallocate(SMAP_AC72_rzmc_struc(n)%rlat)
     if (allocated(SMAP_AC72_rzmc_struc(n)%rlon)) deallocate(SMAP_AC72_rzmc_struc(n)%rlon)
     if (allocated(SMAP_AC72_rzmc_struc(n)%n11 )) deallocate(SMAP_AC72_rzmc_struc(n)%n11 )

     allocate(SMAP_AC72_rzmc_struc(n)%rlat(nlat*nlon))
     allocate(SMAP_AC72_rzmc_struc(n)%rlon(nlat*nlon))
     allocate(SMAP_AC72_rzmc_struc(n)%n11 (nlat*nlon))

     do r=1,nlat
        do c=1,nlon
           SMAP_AC72_rzmc_struc(n)%rlat(c+(r-1)*nlon) = latv(r)
           SMAP_AC72_rzmc_struc(n)%rlon(c+(r-1)*nlon) = lonv(c)
        end do
     end do

     gridDesci = 0.0
     gridDesci(1)  = 0
     gridDesci(2)  = nlon
     gridDesci(3)  = nlat
     gridDesci(4)  = lat0
     gridDesci(5)  = lon0
     gridDesci(6)  = 128
     gridDesci(7)  = latN
     gridDesci(8)  = lonN
     gridDesci(9)  = dlat
     gridDesci(10) = dlon
     gridDesci(20) = 64

     call neighbor_interp_input_withgrid( gridDesci, &
          LIS_rc%obs_gridDesc(k,:), &
          nlon*nlat, &
          SMAP_AC72_rzmc_struc(n)%rlat, &
          SMAP_AC72_rzmc_struc(n)%rlon, &
          SMAP_AC72_rzmc_struc(n)%n11 )
  end if

  !----------------------
  ! Flatten, QC mask (NaN-aware; keep zeros), and interpolate to LIS obs grid
  !----------------------
  allocate(src_data(nlat*nlon))
  allocate(src_mask(nlat*nlon))
  do r=1,nlat
     do c=1,nlon
        src_data(c+(r-1)*nlon) = sm2(r,c)
        is_nan = .not.(sm2(r,c) == sm2(r,c))     ! NaN check
        if ( is_nan .or. ( .not. has_fill_nan .and. sm2(r,c) == fillv ) .or. (sm2(r,c) < 0.0) .or. (sm2(r,c) > 1.0) ) then
           src_mask(c+(r-1)*nlon) = .false.
           src_data(c+(r-1)*nlon) = LIS_rc%udef
        else
           src_mask(c+(r-1)*nlon) = .true.
        end if
     end do
  end do

  ! Quick diagnostics on source field
  n_valid_src = count(src_mask)
  if (n_valid_src > 0) then
     smin = minval(src_data, mask=src_mask)
     smax = maxval(src_data, mask=src_mask)
  else
     smin = 0.0; smax = 0.0
  end if
  write(LIS_logunit,'(A,I0)') 'SMAP src valid pts: ', n_valid_src
  write(LIS_logunit,'(A,2F8.3)') 'src lat span: ', minval(latv), maxval(latv)
  write(LIS_logunit,'(A,2F8.3)') 'src lon span: ', minval(lonv), maxval(lonv)
  write(LIS_logunit,'(A,2(1X,ES13.6))') 'src min/max (valid):', smin, smax

  ! Interp, using a SEPARATE destination mask
  allocate(dst_mask(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
  dst_mask = .false.

  call neighbor_interp( LIS_rc%obs_gridDesc(k,:), &
       src_mask, src_data, &
       dst_mask, smobs_ip, &
       nlon*nlat, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       SMAP_AC72_rzmc_struc(n)%rlat, SMAP_AC72_rzmc_struc(n)%rlon, &
       SMAP_AC72_rzmc_struc(n)%n11, LIS_rc%udef, ios )
  call LIS_verify(ios, 'neighbor_interp failed')

  n_valid_dst = count(dst_mask)
  if (n_valid_dst > 0) then
     dmin = minval(smobs_ip, mask=dst_mask)
     dmax = maxval(smobs_ip, mask=dst_mask)
  else
     dmin = 0.0; dmax = 0.0
  end if
  write(LIS_logunit,'(A,I0)') 'LIS dst valid pts after interp: ', n_valid_dst
  write(LIS_logunit,'(A,2(1X,ES13.6))') 'dst min/max (valid):', dmin, dmax

  ! Cleanup
  deallocate(latv, lonv, sm2, src_data, src_mask, dst_mask)

#else
  ! NetCDF not enabled at compile-time
  smobs_ip = LIS_rc%udef
#endif

contains
  subroutine reverse(v); real, intent(inout) :: v(:); integer i, j; real tmp
    i=1; j=size(v)
    do while(i<j); tmp=v(i); v(i)=v(j); v(j)=tmp; i=i+1; j=j-1; end do
  end subroutine reverse

  subroutine flip_dim1(a); real, intent(inout) :: a(:,:); integer i, j, n1, n2
    n1=size(a,1); n2=size(a,2)
    do j=1,n2
      do i=1,n1/2
        call swap(a(i,j), a(n1-i+1,j))
      end do
    end do
  end subroutine flip_dim1

  subroutine flip_dim2(a); real, intent(inout) :: a(:,:); integer i, j, n1, n2
    n1=size(a,1); n2=size(a,2)
    do i=1,n1
      do j=1,n2/2
        call swap(a(i,j), a(i,n2-j+1))
      end do
    end do
  end subroutine flip_dim2

  subroutine swap(x,y); real, intent(inout) :: x, y; real t; t=x; x=y; y=t; end subroutine swap


end subroutine read_SMAP_AC72_data
!EOP


subroutine create_SMAP_AC72rzmc_filename(ndir, yr, mo, da, filename)
! !USES:
  use LIS_logMod
  implicit none
! !ARGUMENTS:
  character(len=* ) :: ndir
  integer           :: yr, mo, da
  character(len=* ) :: filename
  character (len=4) :: fyr
  character (len=2) :: fmo, fda
  character(len=256) :: cand
  logical :: ex

  write(fyr,'(i4.4)') yr; write(fmo,'(i2.2)') mo; write(fda,'(i2.2)') da
  filename = ''

  cand = trim(ndir)//'/SMAP_RZMC_AC72_'//trim(fyr)//trim(fmo)//trim(fda)//'.nc'
  inquire(file=cand, exist=ex); if (ex) then; filename=cand; return; endif

end subroutine create_SMAP_AC72rzmc_filename

