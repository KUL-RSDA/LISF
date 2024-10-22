!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!
!
! This reader can be used for netcdf files in the following format:
! name: 'SD_yyyymmdd_.nc' with variables 'SD' (m), 'lat' and 'lon'
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_ML_SNWD
! \label{read_ML_SNWD}
!
! !REVISION HISTORY:
!  29 Aug 2019: Hans Lievens; Initial Specification
!
! !INTERFACE: 
subroutine read_ML_SNWD(n,k,OBS_State,OBS_Pert_State) 
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use map_utils
  use LIS_DAobservationsMod
  use LIS_pluginIndices, only : LIS_ML_SNWD_obsId
  use ML_SNWD_Mod, only : ML_SNWD_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  This routine reads Sentinel-1 snow depth observations in netcdf4 format
!  The data is read at 0z every day and is kept in memory. At 
!  each timestep, a subset of data is chosen for use in DA if 
!  the local time of the grid point is 6AM. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[k]    index of the assimilation instance
!  \item[OBS\_State] observations state
!  \item[OBS\_Pert\_State] observations perturbation state
!  \end{description}
!
!EOP
  type(ESMF_Field)              :: snowField,pertfield
  logical                       :: alarmCheck
  logical                       :: data_upd, file_exists
  logical                       :: dataflag(LIS_npes)
  logical                       :: dataflag_local
  integer                       :: c,r, p, t
  character*100                 :: obsdir
  character*80                  :: ML_filename
  integer                       :: ftn
  real                          :: lon, lhour
  real                          :: gmt
  real                          :: dt
  integer                       :: zone
  integer                       :: grid_index
  real                          :: ssdev(LIS_rc%obs_ngrid(k))
  real,             pointer     :: obsl(:)
  integer                       :: gid(LIS_rc%obs_ngrid(k))
  integer                       :: assimflag(LIS_rc%obs_ngrid(k))
  integer                       :: status, iret, ierr
  integer                       :: fnd
  real                          :: snwd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                       :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/  !BZ



  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status,'Error in AttributeGet: Data Directory')

!-------------------------------------------------------------------------
!   Read the data at 0z daily. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ML snow depth read alarm")

  if(alarmCheck.or.ML_SNWD_struc(n)%startMode) then 
     ML_SNWD_struc(n)%startMode = .false.
     
     ML_SNWD_struc(n)%snwd = LIS_rc%udef
     ML_SNWD_struc(n)%snwdtime = -1

     call ML_SNWD_filename(ML_filename,obsdir,&
          LIS_rc%yr,LIS_rc%mo,LIS_rc%da)       

     inquire(file=ML_filename,exist=file_exists)
     if(file_exists) then 

        write(LIS_logunit,*)  '[INFO] Reading ML SNWD data ',ML_filename
        call read_ML_SNWD_data(n,k, ML_filename, ML_SNWD_struc(n)%snwd)


!-------------------------------------------------------------------------
! Store the GMT corresponding to 6AM localtime at each grid point
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
                 lon = LIS_obs_domain(n,k)%lon(grid_index)
                 
                 lhour = 6.0
                 call LIS_localtime2gmt(gmt,lon,lhour,zone)
                 ML_SNWD_struc(n)%snwdtime(c,r) = gmt

              endif
           enddo
        enddo

     endif
  endif

!-------------------------------------------------------------------------
! Update the OBS_State
!-------------------------------------------------------------------------

  call ESMF_StateGet(OBS_State,"Observation01",snowfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(snowfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  obsl = LIS_rc%udef 

!-------------------------------------------------------------------------
!  Update the OBS_State by subsetting to the local grid time  
!-------------------------------------------------------------------------     

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)
           
           dt = (LIS_rc%gmt - ML_SNWD_struc(n)%snwdtime(c,r))*3600.0
           lon = LIS_obs_domain(n,k)%lon(grid_index)

           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = & 
                   ML_SNWD_struc(n)%snwd(c,r)
           endif
           
        endif
     enddo
  enddo

  dataflag_local = .false. 

!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------     

  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_ML_SNWD_obsId)//char(0), & 
       n, k,OBS_state)

  snwd_current = LIS_rc%udef
  call LIS_checkForValidObs(n,k,obsl,fnd,snwd_current)

  if(fnd.eq.0) then 
     dataflag_local = .false. 
  else
     dataflag_local = .true. 
  endif
 
#if (defined SPMD)
  call MPI_ALLGATHER(dataflag_local,1, MPI_LOGICAL, dataflag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, ierr)
#endif
  data_upd = .false.
  
  do p=1,LIS_npes
     data_upd = data_upd.or.dataflag(p)
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
          .true., rc=status)
     call LIS_verify(status, 'Error: AttributeSet in Data Update Status')
     

     call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
          rc=status)
     call LIS_verify(status, 'ESMF_StateGet for Observation01 for OBS_Pert_State failed in read_ML_SNWD')
     
     if(LIS_rc%obs_ngrid(k).gt.0) then 

!linearly scale the observation err
        ssdev = ML_SNWD_struc(n)%ssdev 
        do t=1,LIS_rc%obs_ngrid(k)
           if(obsl(t).ne.-9999.0) then 
              ssdev(t) =  ML_SNWD_struc(n)%ssdev !+ 0.05*obsl(t)
           endif
        enddo

        call ESMF_AttributeSet(pertField,"Standard Deviation",&
             ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

! To be added (for EnKS test):
!        call ESMF_AttributeSet(OBS_State,&
!             name="Data averaging factor",&
!             value=float(days(LIS_rc%mo)),rc=status)
!        call LIS_verify(status)

        call ESMF_AttributeSet(snowfield,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status,'Error: AttributeSet in Grid Number')
        
        call ESMF_AttributeSet(snowfield,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status, 'Error: AttributeSet in Assimilation Flag')
        
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status, "Error: AttributeSet Data Update Status")
     return
  end if
 

end subroutine read_ML_SNWD



!BOP
!
! !ROUTINE: read_ML_SNWD_data
! \label{read_ML_SNWD_data}
!
! !INTERFACE:
subroutine read_ML_SNWD_data(n, k, fname, snwd_ip)
!
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod
  use map_utils,    only : latlon_to_ij
  use ML_SNWD_Mod, only : ML_SNWD_struc

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
!  This subroutine reads the ML NETCDF files
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ML SNWD file
!  \item[SDobs\_ip]   snow depth data processed to the LIS domain
! \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  real                        :: SD(ML_SNWD_struc(n)%nr,ML_SNWD_struc(n)%nc)
  real                        :: lat_nc(ML_SNWD_struc(n)%nr)
  real                        :: lon_nc(ML_SNWD_struc(n)%nc)
  real                        :: snwd_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                     :: nsnwd_ip(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  logical                     :: file_exists
  integer                     :: c,r,i,j
  integer                     :: stn_col,stn_row
  real                        :: col,row
  integer                     :: nid
  integer                     :: SDId,latId,lonId
  integer                     :: ios

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

!  SD = LIS_rc%udef
!  lat_nc = LIS_rc%udef
!  lon_nc = LIS_rc%udef
!  snwd_ip = LIS_rc%udef

  inquire(file=fname, exist=file_exists)
  if(file_exists) then
     write(LIS_logunit,*) 'Reading ',trim(fname)
     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error opening file '//trim(fname))
 
     ! variables
     ios = nf90_inq_varid(nid, 'SD',SDid)
     call LIS_verify(ios, 'Error nf90_inq_varid: snow depth data')

     ios = nf90_inq_varid(nid, 'lat',latid)
     call LIS_verify(ios, 'Error nf90_inq_varid: latitude data')

     ios = nf90_inq_varid(nid, 'lon',lonid)
     call LIS_verify(ios, 'Error nf90_inq_varid: longitude data')

     !values
     ios = nf90_get_var(nid, SDid, SD)
     call LIS_verify(ios, 'Error nf90_get_var: SD')

     ios = nf90_get_var(nid, latid, lat_nc)
     call LIS_verify(ios, 'Error nf90_get_var: lat')

     ios = nf90_get_var(nid, lonid, lon_nc)
     call LIS_verify(ios, 'Error nf90_get_var: lon')

     ! close file
     ios = nf90_close(ncid=nid)
     call LIS_verify(ios,'Error closing file '//trim(fname))



     snwd_ip = 0
     nsnwd_ip = 0

     ! Map the data by averaging 
     do i=1,ML_SNWD_struc(n)%nr
        do j=1,ML_SNWD_struc(n)%nc
              
           call latlon_to_ij(LIS_domain(n)%lisproj,&
                lat_nc(i),lon_nc(j),col,row)
           stn_col = nint(col)
           stn_row = nint(row)

           if(SD(i,j).ge.-999.and.&
                stn_col.gt.0.and.stn_col.le.LIS_rc%obs_lnc(k).and.&
                stn_row.gt.0.and.stn_row.le.LIS_rc%obs_lnr(k)) then
              snwd_ip(stn_col,stn_row) = snwd_ip(stn_col,stn_row) + SD(i,j)
              nsnwd_ip(stn_col,stn_row) = nsnwd_ip(stn_col,stn_row) + 1
           endif
        enddo
     enddo

     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(nsnwd_ip(c,r).ne.0) then
              !Sentinel-1 SD bias correction factor 1.0/0.85 if wanted 
              ! is implemented here, now disabled
              snwd_ip(c,r) = 1.0*snwd_ip(c,r)/nsnwd_ip(c,r) ! NoahMP4 snow DA wants it in meter
           else
              snwd_ip(c,r) = LIS_rc%udef
           endif
        enddo
     enddo


  endif

#endif

end subroutine read_ML_SNWD_data



!BOP
!
! !ROUTINE: ML_SNWD_filename
! \label{ML_SWND_filename}
! 
! !INTERFACE: 
subroutine ML_SNWD_filename(filename, ndir, yr, mo, da)
  
  implicit none
! !ARGUMENTS: 
  character*80      :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped ML filename
!  
!  The arguments are: 
!  \begin{description}
!  \item[name] name of the ML filename
!  \item[ndir] name of the ML root directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da 

  filename = trim(ndir)//'/SD_'//trim(fyr)//trim(fmo)//trim(fda)//'_.nc'
    
end subroutine ML_SNWD_filename



