!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_SMOSL2sm
! \label{read_SMOSL2sm}
!
! !REVISION HISTORY:
!  16 Dec 14    Sujay Kumar; Initial specification
!  09 Jan 24    Sara Modanesi; Updated subroutines to activate the DA component and ...
!               read ncfiles instead of .DBL format
! !INTERFACE: 

!check SM09012024 added the k in the function
subroutine read_SMOSL2sm(n, k, OBS_State, OBS_Pert_State)
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
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use SMOSL2sm_Mod, only : SMOSL2sm_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k !SM09012024
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real,  parameter       :: MAX_SM_VALUE=0.45, MIN_SM_VALUE=0.0001
  integer                :: ftn
  integer                :: ierr,status
  integer                :: grid_index
  character(len=LIS_CONST_PATH_LEN) :: smobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname_A, fname_D
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p
  real,          pointer :: obsl(:)
  type(ESMF_Time)        :: obsStartTime,obsEndTime,obsTime
  type(ESMF_Time)        :: currTime, prevTime
  type(ESMF_TimeInterval) :: dayInt,dt,lis_dt, zero_dt
  type(ESMF_Field)       :: smfield
  integer                :: sind
  character*100          :: temp1
  character*1            :: fproc(4)
  character(len=LIS_CONST_PATH_LEN) :: smos_filename
  integer                :: yr,mo,da,hr,mn,ss
  character*8            :: yyyymmdd
  character(len=LIS_CONST_PATH_LEN) :: list_files
  ! SM09012024 added LIS_rc%obs_nigrid(k) instead LIS_rc%ngrid(n) and modified with LIS_rc%obs_lnc(k) and LIS_rc%obs_lnr(k)
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: smobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: sm_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real                   :: smvalue
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
  type(ESMF_Time)        :: acTime

  integer                :: ts_check,yr_t,mo_t,da_t, hr_t, mn_t,ss_t

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       smobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call ESMF_TimeintervalSet(dayInt, d=1,rc=status)
  call LIS_verify(status, 'ESMF_TimeIntervalSet failed in read_SMOSL2sm')

  data_upd = .false. 
  obs_unsc = LIS_rc%udef
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "SMOS L2 read alarm")
  
  if(alarmCheck.or.SMOSL2sm_struc(n)%startMode) then 
     SMOSL2sm_struc(n)%startMode = .false.

     SMOSL2sm_struc(n)%smobs = LIS_rc%udef
     call ESMF_TimeSet(acTime, yy=2000,&
          mm = 1, dd=1, h =0, &
          m = 0, s = 0, calendar = LIS_calendar,&
          rc=status)
     call LIS_verify(status, 'ESMF_TimeSet failed in read_SMOSL2_data (3)')
     SMOSL2sm_struc(n)%smTime = acTime

     write(unit=temp1,fmt='(i4.4)') LIS_localPet
     read(unit=temp1,fmt='(4a1)')fproc
     
! dump the list of files for the current date to a file (note that
! we assume a flat organization of the files under the SMOS observation
! directory. 
     !SM19012024
     write(yyyymmdd,'(i4.4,2i2.2)') LIS_rc%yr, LIS_rc%mo, LIS_rc%da
     list_files = 'ls '//trim(SMOSL2sm_struc(n)%odir)//'/*'//trim(yyyymmdd)&
          //"*.nc > SMOS_filelist."//&
          fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat"
     call system(trim(list_files))

     ftn = LIS_getNextUnitNumber()
     open(ftn,file="./SMOS_filelist."//&
          fproc(1)//fproc(2)//fproc(3)//fproc(4)//".dat",status='old',iostat=ierr)

     do while(ierr.eq.0) 
        read(ftn,'(a)',iostat=ierr) smos_filename
        if(ierr.ne.0) then 
           exit
        else
! check first if the observation time is within the LIS time window (daily)
           call ESMF_TimeSet(currTime, yy=LIS_rc%yr,&
                mm = LIS_rc%mo, dd=LIS_rc%da, h = LIS_rc%hr, &
                m = LIS_rc%mn, s = LIS_rc%ss, calendar = LIS_calendar,&
                rc=status)
           call LIS_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (1)')
           prevTime = currTime 
           currTime = currTime + dayInt
           
           sind = index(smos_filename, "SMUDP2_")
           sind = sind+7
           read(smos_filename(sind:sind+3),'(i4)') yr
           read(smos_filename(sind+4:sind+5),'(i2)') mo
           read(smos_filename(sind+6:sind+7),'(i2)') da
           read(smos_filename(sind+9:sind+10),'(i2)') hr
           read(smos_filename(sind+11:sind+12),'(i2)') mn
           read(smos_filename(sind+13:sind+14),'(i2)') ss
           
           call ESMF_TimeSet(obsStartTime, yy=yr,&
                mm = mo, dd=da, h = hr, &
                m = mn, s = ss, calendar = LIS_calendar,&
                rc=status)
           call LIS_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (2)')

           read(smos_filename(sind+16:sind+19),'(i4)') yr
           read(smos_filename(sind+20:sind+21),'(i2)') mo
           read(smos_filename(sind+22:sind+23),'(i2)') da
           read(smos_filename(sind+25:sind+26),'(i2)') hr
           read(smos_filename(sind+27:sind+28),'(i2)') mn
           read(smos_filename(sind+29:sind+30),'(i2)') ss

           call ESMF_TimeSet(obsEndTime, yy=yr,&
                mm = mo, dd=da, h = hr, &
                m = mn, s = ss, calendar = LIS_calendar,&
                rc=status)
           call LIS_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (3)')
 
           dt = (obsEndTime - obsStartTime)/2
           obsTime = (obsStartTime + dt)

           if((obsTime.gt.prevTime).and.(obsTime.le.currTime)) then 
              call read_SMOSL2_data(n,k, smos_filename)
           endif
        endif
     enddo

     call LIS_releaseUnitNumber(ftn)

  endif



  fnd = 0 
  sm_current = LIS_rc%udef

  ! dt is not defined as absolute value of the time difference to avoid
  ! double counting of the data in assimilation. 

  call ESMF_TimeSet(currTime, yy=LIS_rc%yr,&
       mm = LIS_rc%mo, dd=LIS_rc%da, h = LIS_rc%hr, &
       m = LIS_rc%mn, s = LIS_rc%ss, calendar = LIS_calendar,&
       rc=status)
  call LIS_verify(status, 'ESMF_TimeSet failed in readSMOSL2smObs (1)')
  call ESMF_TimeIntervalSet(lis_dt, s=nint(LIS_rc%ts))
  call ESMF_TimeIntervalSet(zero_dt, s=0)
  !SM09012024 changed LIS_rc%lnr(n) with LIS_rc%obs_lnr(k). Same for LIS_rc%lnc(n) 
  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           dt = (currTime - SMOSL2sm_struc(n)%smtime(c,r))
           if(dt.ge.zero_dt.and.dt.lt.lis_dt) then 
              sm_current(c,r) = & 
                   SMOSL2sm_struc(n)%smobs(c,r)
              fnd = 1
           endif
        endif
     enddo
  enddo



  !-------------------------------------------------------------------------
  !  Transform data to the LSM climatology using a CDF-scaling approach
  !-------------------------------------------------------------------------     
  !SM03052024 This was modified to exclude the storing of unscaled obs from the if
  !i.e. the store is done even if the rescaling is not called
  do r =1,LIS_rc%obs_lnr(k)
     do c =1,LIS_rc%obs_lnc(k)
        if (LIS_domain(n)%gindex(c,r) .ne. -1)then
           obs_unsc(LIS_domain(n)%gindex(c,r)) = &
                sm_current(c,r)
        end if
     end do
  end do

  if(SMOSL2sm_struc(n)%scal.ne.0.and.fnd.ne.0) then
     !SM09012024 added k in the function
     call LIS_rescale_with_CDF_matching(    &
          n,k,                                   & 
          SMOSL2sm_struc(n)%nbins,         & 
          SMOSL2sm_struc(n)%ntimes,        & 
          MAX_SM_VALUE,                        & 
          MIN_SM_VALUE,                        & 
          SMOSL2sm_struc(n)%model_xrange,  &
          SMOSL2sm_struc(n)%obs_xrange,    &
          SMOSL2sm_struc(n)%model_cdf,     &
          SMOSL2sm_struc(n)%obs_cdf,       &
          sm_current)
  endif

!-------------------------------------------------------------------------
!  LSM quality control on obs 
!-------------------------------------------------------------------------     

  call ESMF_StateGet(OBS_State,"Observation01",smfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')

  call ESMF_FieldGet(smfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  !SM09012024
  obsl = LIS_rc%udef 
  if(fnd.ne.0) then 
     do r=1, LIS_rc%obs_lnr(k)
        do c=1, LIS_rc%obs_lnc(k)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              obsl(LIS_domain(n)%gindex(c,r))=sm_current(c,r)
           endif
        enddo
     enddo
  endif

   !-------------------------------------------------------------------------
   !  Apply LSM based QC and screening of observations
   !-------------------------------------------------------------------------
  !SM05052024 modified based on the NASA_SMAPsm reader
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_SMOSL2smobsId)//char(0),n,k, OBS_state) 

  call LIS_checkForValidObs(n, k, obsl, fnd, sm_current) 

   if (fnd .eq. 0) then
      data_upd_flag_local = .false.
   else
      data_upd_flag_local = .true.
   endif

#if (defined SPMD)
   call MPI_ALLGATHER(data_upd_flag_local, 1, &
                      MPI_LOGICAL, data_upd_flag(:), &
                      1, MPI_LOGICAL, LIS_mpi_comm, status)
   data_upd = any(data_upd_flag)
#else
   data_upd = data_upd_flag_local
#endif


  !SM09012024 all the LIS_rc%ngrid(n) where changed to LIS_rc%obs_ngrid(k)
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
        call ESMF_AttributeSet(smField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(smField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(smfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_SMOSL2sm

!BOP
! 
! !ROUTINE: read_SMOSL2_data
! \label{read_SMOSL2_data}
!
! !INTERFACE:

!SM01092024 'k' added + change of LIS_rc%obs_lnc(n) and LIS_rc%obs_lnr(n)
!SM01092024 changed to read nc files + updated flags of SMOS obs (check confidence and science flags)
subroutine read_SMOSL2_data(n, k, fname)
! 
! !USES:  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use map_utils
  use SMOSL2sm_Mod, only : SMOSL2sm_struc

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n
  integer                       :: k 
  character (len=*)             :: fname

! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads a given SMOS L2 soil moisture data, applies the
!  quality control flags and then grids the data into the LIS/LIS grid. 
!   
!  The data format is described in the SMOS Level 2 auxilary data products
!  specifications document (pages 78 onwards)
!  \hyperref{https://earth.esa.int/documents/10174/127856/SO-TN-IDR-GS-0006_L2+Spec_v6.1_2012-02-09.pdf/ea56dca3-fe04-4196-a200-ab34da105a9c?version=1.0}{}{}{https://earth.esa.int/documents/10174/127856/SO-TN-IDR-GS-0006\_L2+Spec\_v6.1\_2012-02-09.pdf/ea56dca3-fe04-4196-a200-ab34da105a9c?version=1.0}
!  
!  The code here applies available dataflags in the data for inline
!  quality control. Data is excluded when the following criteria is met
!  \begin{itemize}
!   \item RFI probabilty is $>$ 50 
!   \item vegetation optical thickness is $>$ 0.8
!   \item data quality is reported to be poor
!   \item when reported uncertainty is high
!   SM01092023 
!   \item considering confidence and science flags
!  \end{itemize}
!
!  For grid points with multiple data values during a day, data at the
!  latest observation time is chosen. The observation times in the 
!  SMOS data is reported as the number of days since Jan 1, 2000. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOSL2 file
!  \item[smobs]    soil moisture data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
  !integer                          :: ftn
  integer                          :: nid
  integer                          :: sm_ID,lat_ID,lon_ID,sm_dqx_ID,opt_thick_ID,tsurf_ID !totpts_ID
  integer                          :: confidence_flags_ID, science_flags_ID
  integer                          :: gqx_ID, grid_point_id_ID, acday_ID, acsec_ID, acms_ID, rfi_prob_ID
  real, pointer                    :: sm(:), lat(:), lon(:), sm_dqx(:), opt_thick(:), tsurf(:), rfi_prob(:), gqx(:)
  integer, pointer                 :: confidence_flags(:), science_flags(:)
  integer                          :: lon_size, dimid
  integer                          :: i,c,r
  real                             :: col,row
  integer, pointer                 :: grid_point_id(:)
  integer, pointer                 :: acday(:),acsec(:),acms(:)
  integer                          :: yr,mo,da,da1,hr,mn,ss
  logical                          :: file_exists
  integer                          :: stc, enc, str, enr, c1,r1
  integer                          :: status
  type(ESMF_Time)                  :: acTime
  integer                          :: ios
  logical                          :: cond1, cond2, cond3, cond4

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=fname,exist=file_exists)
  
  if(file_exists) then 
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
     
     ! variables (in_varid + dimensions)
     ios = nf90_inq_varid(nid, 'Longitude',lon_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Longitude')

     ios=nf90_inq_dimid(nid,'n_grid_points',dimid)
     call LIS_verify(ios, 'nf90_inq_dimid: Longitude')

     ios = nf90_inquire_dimension(nid, dimid, len=lon_size)
     call LIS_verify(ios, 'nf90_inquire_dimension: Longitude')
     
     !-latitude-!
     ios = nf90_inq_varid(nid, 'Latitude',lat_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Latitude')

     !-soil moisture-!
     ios = nf90_inq_varid(nid, 'Soil_Moisture',sm_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Soil_Moisture')

     !-soil moisture dqx-! 
     ios = nf90_inq_varid(nid, 'Soil_Moisture_DQX',sm_dqx_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Soil_Moisture_DQX')

     !-optical thickness-!
     ios = nf90_inq_varid(nid, 'Optical_Thickness_Nad',opt_thick_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Optical_Thickness_Nad')

     !-surface temperature-!
     ios = nf90_inq_varid(nid, 'Surface_Temperature',tsurf_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Surface_Temperature')

     !-GQX-!
     ios = nf90_inq_varid(nid, 'GQX',gqx_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: GQX')

     !-grid points ID-!
     ios = nf90_inq_varid(nid, 'Grid_Point_ID',grid_point_id_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Grid_Point_ID')

     !-Days-!
     ios = nf90_inq_varid(nid, 'Days',acday_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Days')

     !-seconds-!
     ios = nf90_inq_varid(nid, 'Seconds',acsec_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Seconds')

     !-Microseconds-!
     ios = nf90_inq_varid(nid, 'Microseconds',acms_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Microseconds')

     !-RFI probability-!
     ios = nf90_inq_varid(nid, 'RFI_Prob',rfi_prob_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: RFI_Prob')

     !-Confidence Flags-!
     ios = nf90_inq_varid(nid, 'Confidence_Flags', confidence_flags_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Confidence_Flags')

     !-Science Flags-!
     ios=nf90_inq_varid(nid, 'Science_Flags', science_flags_ID)
     call LIS_verify(ios, 'Error nf90_inq_varid: Science_Flags')

     allocate(lat(lon_size))
     ios = nf90_get_var(nid, lat_ID, lat)
     call LIS_verify(ios, 'Error nf90_get_var: Latitude')

     allocate(lon(lon_size))
     ios = nf90_get_var(nid, lon_ID, lon)
     call LIS_verify(ios, 'Error nf90_get_var: Longitude')

     allocate(sm(lon_size))
     ios = nf90_get_var(nid, sm_ID, sm)
     call LIS_verify(ios, 'Error nf90_get_var: Soil_Moisture')

     allocate(sm_dqx(lon_size))
     ios = nf90_get_var(nid, sm_dqx_ID, sm_dqx)
     call LIS_verify(ios, 'Error nf90_get_var: Soil_Moisture_DQX')

     allocate(opt_thick(lon_size))
     ios = nf90_get_var(nid, opt_thick_ID, opt_thick)
     call LIS_verify(ios, 'Error nf90_get_var: Optical_Thickness_Nad')

     allocate(tsurf(lon_size))
     ios = nf90_get_var(nid, tsurf_ID, tsurf)
     call LIS_verify(ios, 'Error nf90_get_var: Surface_Temperature')

     allocate(gqx(lon_size))
     ios = nf90_get_var(nid, gqx_ID, gqx)
     call LIS_verify(ios, 'Error nf90_get_var: GQX')

     allocate(grid_point_id(lon_size))
     ios = nf90_get_var(nid, grid_point_id_ID, grid_point_id)
     call LIS_verify(ios, 'Error nf90_get_var: Grid_Point_ID')
     
     allocate(acday(lon_size))
     ios = nf90_get_var(nid, acday_ID, acday)
     call LIS_verify(ios, 'Error nf90_get_var: Days')

     allocate(acsec(lon_size))
     ios = nf90_get_var(nid, acsec_ID, acsec)
     call LIS_verify(ios, 'Error nf90_get_var: Seconds')

     allocate(acms(lon_size))
     ios = nf90_get_var(nid, acms_ID, acms)
     call LIS_verify(ios, 'Error nf90_get_var: Microseconds')

     allocate(rfi_prob(lon_size))
     ios = nf90_get_var(nid, rfi_prob_ID, rfi_prob)
     call LIS_verify(ios, 'Error nf90_get_var: RFI_prob')

     allocate(confidence_flags(lon_size))
     ios = nf90_get_var(nid, confidence_flags_ID, confidence_flags)
     call LIS_verify(ios, 'Error nf90_get_var: Confidence_Flags')

     allocate(science_flags(lon_size))
     ios = nf90_get_var(nid, science_flags_ID, science_flags)
     call LIS_verify(ios, 'Error nf90_get_var: Science_Flags')

     ! close file
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))
  
     do i=1, size(lon,1) !totpts

        if((.not.isNaN(lat(i))).and.(.not.isNaN(lon(i))).and.&
             abs(lat(i)).lt.400.and.abs(lon(i)).lt.400) then
          call latlon_to_ij(LIS_domain(n)%lisproj,lat(i),lon(i),&
               col,row)

          c = nint(col)
          r = nint(row)
          stc = max(1,c)
          enc = min(LIS_rc%obs_lnc(k),c)
          str = max(1,r)
          enr = min(LIS_rc%obs_lnr(k),r)
          !SM04032024 added conditions for science and confidence flags          
          cond1=(sm(i).gt.0.001.and.sm(i).lt.1.00).and.(sm_dqx(i).le.0.1).and.(gqx(i).lt.10).and.(rfi_prob(i).lt.50).and.(opt_thick(i).lt.0.8)
          cond2=(.not. btest(confidence_flags(i), 1)).and.(.not. btest(confidence_flags(i), 2)).and.&
                  (.not. btest(confidence_flags(i), 4)).and.(.not. btest(confidence_flags(i), 5)).and.(.not. btest(confidence_flags(i), 6))
          cond3=(.not. btest(science_flags(i), 0)).and.(.not. btest(science_flags(i), 5)).and.(.not. btest(science_flags(i), 16)).and.&
                  (.not. btest(science_flags(i),18)).and.(.not. btest(science_flags(i), 19)).and.&
                  (.not. btest(science_flags(i), 26))

          do c1=stc,enc
             do r1=str,enr
             cond4=(c1.ge.1.and.c1.le.LIS_rc%obs_lnc(k)).and.(r1.ge.1.and.r1.le.LIS_rc%obs_lnr(k))
                if (cond1 .and. cond2 .and. cond3 .and. cond4) then

                   if(SMOSL2sm_struc(n)%smobs(c1,r1).gt.0) then
                     !              print*, 'data already there',c,r
                   else
                      if(SMOSL2sm_struc(n)%smobs(c1,r1).eq.-9999.0) then !if data is not present already
                         call SMOS_julhr_date( hr, da, mo, yr, acday(i)*24 ) 
                         call LIS_seconds2time(acsec(i),da1, hr, mn, ss)
                         call ESMF_TimeSet(acTime, yy=yr,&
                              mm = mo, dd=da, h = hr, &
                              m = mn, s = ss, calendar = LIS_calendar,&
                              rc=status)
                         call LIS_verify(status, 'ESMF_TimeSet failed in read_SMOSL2_data (4)')
                         if(acTime > SMOSL2sm_struc(n)%smTime(c1,r1)) then 
                            SMOSL2sm_struc(n)%smTime(c1,r1) = acTime                       
                            SMOSL2sm_struc(n)%smobs(c1,r1) = sm (i)
                         endif
                      endif
                   endif
                endif 
             enddo
          enddo
        endif
     enddo 
  endif
  deallocate(lat)
  deallocate(lon)
  deallocate(sm)
  deallocate(sm_dqx)
  deallocate(opt_thick)
  deallocate(tsurf)
  deallocate(gqx)
  deallocate(grid_point_id)
  deallocate(acday)
  deallocate(acsec)
  deallocate(acms)
  deallocate(rfi_prob)
  deallocate(confidence_flags)
  deallocate(science_flags)
#endif

end subroutine read_SMOSL2_data

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! 
! !INTERFACE:    
subroutine SMOS_julhr_date( hour, day, month, year, julhr ) 

  implicit none 
! !ARGUMENTS:   
  integer,  intent(out)       :: day  
  integer,  intent(out)       :: hour 
  integer,  intent(out)       :: month
  integer,  intent(out)       :: year 
  integer,  intent(in)        :: julhr
!
! !DESCRIPTION:
!     uses the julian hour to determine the 2-digit zulu time, day,  
!     month, and 4-digit year.    
!
!     
!    \subsubsection{Method}
!     - determine the current zulu hour   \newline
!     - determine the total number of elapsed days   \newline
!     - count forward to the current day/month/year   \newline
!
!    The arguments are: 
!    \begin{description}
!    \item[hour]
!     the zulu time of the julian hour   
!    \item[day]
!     day of the month (1..31)   
!    \item[month]
!     month of the year (1..12)  
!    \item[year]
!     four digit year
!    \item[julhr]
!     the julian hour being processed
!    \end{description}
!
!EOP  
  integer                     :: dypmon(12)   
  integer                     :: elapdy   
  integer                     :: i
  logical                     :: done 
  
  data dypmon   /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  
!     ------------------------------------------------------------------
!     initialize done flag to false.
!     ------------------------------------------------------------------

  done = .false.
      
!    ------------------------------------------------------------------
!     extract the zulu hour from the julian hour.   
!     ------------------------------------------------------------------
          
  hour = mod(julhr, 24) 

!    ------------------------------------------------------------------
!     determine the number of days that have elapsed since dec 31, 1967 
!     (julian hour 0).  
!     ------------------------------------------------------------------
    
  elapdy = julhr / 24   
          
!     ------------------------------------------------------------------
!     initialize the starting day, month, and year values.  
!     ------------------------------------------------------------------
    
  if (elapdy .gt. 0) then        
     year   = 2000
     day    = 2
     month  = 1
     elapdy = elapdy - 1      
  else       
     print*, 'error in SMOS_julhr_date'
     stop
  endif
      
!    ------------------------------------------------------------------
!     loop through the elapsed days to determine the year.  
!     ------------------------------------------------------------------
    
  do while(.not.done) 

     dypmon(2) = 28  

!    ------------------------------------------------------------------
!     determine if year value is a leap year.  leap years occur every 
!     4 years, with the exception of century years not evenly 
!     divisible by 400.   
!     ------------------------------------------------------------------
    
     if ((mod(year, 4)   .eq. 0) .and.   &
          ((mod(year, 100) .ne. 0) .or. &
          (mod(year, 400) .eq. 0))) then  
        
        dypmon(2) = 29
        
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is more than a year's worth,  
!     subtract the appropriate number of days, and increment the year 
!     value.  
!     ------------------------------------------------------------------
    
     if (dypmon(2) .eq. 28) then      
        if (elapdy .ge. 365) then         
           year = year + 1 
           elapdy = elapdy - 365           
        else          
           done = .true.   
        endif
     else     
        if (elapdy .ge. 366) then         
           year = year + 1 
           elapdy = elapdy - 366           
        else          
           done = .true.           
        endif
     endif

!     ------------------------------------------------------------------
!     if the elapsed number of days is less than a year's worth, then   
!     exit loop.
!     ------------------------------------------------------------------    
  enddo

!     ------------------------------------------------------------------
!     count the days and months elapsed in the current year.
!     ------------------------------------------------------------------
    
  do i = 1, elapdy      
     day = day + 1
     if (day .gt. dypmon(month)) then        
        day = 1   
        month = month + 1                
        if(month.gt.12) then 
           month = 1
           year = year + 1
        endif
     endif     
  enddo

  return
  
end subroutine SMOS_julhr_date

