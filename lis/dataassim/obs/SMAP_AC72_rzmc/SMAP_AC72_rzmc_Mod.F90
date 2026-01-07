!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: SMAP_AC72_rzmc_Mod
!
! !DESCRIPTION:
! Observation plugin for SMAP AC72 root-zone soil moisture.
! Modified (2025-10-10) to:
!  - stop hard-coding a 0.25° source grid
!  - let the reader discover lat/lon and build the source grid mapping
!
module SMAP_AC72_rzmc_Mod
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  private

  public :: SMAP_AC72_rzmc_setup
  public :: SMAP_AC72_rzmc_struc

  type, public :: SMAP_AC72_rzmc_dec
     logical                :: startMode
     character(len=8)       :: sensor
     real                   :: version
     integer                :: useSsdevScal
     integer                :: nc
     integer                :: nr
     real,     allocatable  :: smobs(:,:)
     real,     allocatable  :: smtime(:,:)

     real                   :: ssdev_inp
     integer                :: ecvnc, ecvnr
     type(proj_info)        :: ecvproj
     integer, allocatable   :: n11(:)
     real,    allocatable   :: rlat(:)
     real,    allocatable   :: rlon(:)

     real,    allocatable   :: model_xrange(:,:,:)
     real,    allocatable   :: obs_xrange(:,:,:)
     real,    allocatable   :: model_cdf(:,:,:)
     real,    allocatable   :: obs_cdf(:,:,:)
     real,    allocatable   :: model_mu(:,:)
     real,    allocatable   :: obs_mu(:,:)
     real,    allocatable   :: model_sigma(:,:)
     real,    allocatable   :: obs_sigma(:,:)

     integer                :: nbins
     integer                :: ntimes
     logical                :: cdf_read_mon
     integer                :: cdf_read_opt
     character*100          :: modelcdffile
     character*100          :: obscdffile
     logical                :: midnight_assimilation
     logical                :: barren_assimilation
  end type SMAP_AC72_rzmc_dec

  type(SMAP_AC72_rzmc_dec), allocatable :: SMAP_AC72_rzmc_struc(:)

contains

!BOP
! !ROUTINE: SMAP_AC72_rzmc_setup
! \label{SMAP_AC72_rzmc_setup}
subroutine SMAP_AC72_rzmc_setup(k, OBS_State, OBS_Pert_State)
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_historyMod
  use LIS_dataAssimMod
  use LIS_perturbMod
  use LIS_DAobservationsMod
  use LIS_logmod
  implicit none
  integer                ::  k
  type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
  type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)

  real, parameter        ::  minssdev = 0.001
  integer                ::  n,i,t,jj,status,ftn,ngrid
  type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
  type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
  type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
  type(ESMF_ArraySpec)   ::  pertArrSpec
  character(len=LIS_CONST_PATH_LEN) ::  esaccismobsdir
  character*100          ::  temp
  character*1            ::  vid(2)
  character*40, allocatable  ::  vname(:)
  real        , allocatable  ::  varmin(:), varmax(:)
  type(pert_dec_type)    ::  obs_pert
  real, pointer          ::  obs_temp(:,:)
  real,      allocatable ::  ssdev(:)

  allocate(SMAP_AC72_rzmc_struc(LIS_rc%nnest))

  call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4, rc=status); call LIS_verify(status)
  call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4, rc=status); call LIS_verify(status)
  call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4, rc=status); call LIS_verify(status)

  call ESMF_ConfigFindLabel(LIS_config,"SMAP AC72 root zone soil moisture data directory:", rc=status)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,esaccismobsdir, rc=status)
     call LIS_verify(status, 'SMAP AC72 root zone soil moisture data directory: is missing')
     call ESMF_AttributeSet(OBS_State(n),"Data Directory", esaccismobsdir, rc=status); call LIS_verify(status)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SMAP AC72 root zone use scaled standard deviation model:", rc=status)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,SMAP_AC72_rzmc_struc(n)%useSsdevScal, rc=status)
     call LIS_verify(status, "SMAP AC72 root zone use scaled standard deviation model: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SMAP AC72 root zone model CDF file:", rc=status)
  do n=1,LIS_rc%nnest
     if (LIS_rc%dascaloption(k).ne."none") then
        call ESMF_ConfigGetAttribute(LIS_config,SMAP_AC72_rzmc_struc(n)%modelcdffile,rc=status)
        call LIS_verify(status, 'SMAP AC72 root zone model CDF file: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SMAP AC72 root zone observation CDF file:", rc=status)
  do n=1,LIS_rc%nnest
     if (LIS_rc%dascaloption(k).ne."none") then
        call ESMF_ConfigGetAttribute(LIS_config,SMAP_AC72_rzmc_struc(n)%obscdffile,rc=status)
        call LIS_verify(status, 'SMAP AC72 root zone observation CDF file: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config, "SMAP AC72 root zone soil moisture number of bins in the CDF:", rc=status)
  do n=1,LIS_rc%nnest
     if (LIS_rc%dascaloption(k).ne."none") then
        call ESMF_ConfigGetAttribute(LIS_config,SMAP_AC72_rzmc_struc(n)%nbins, rc=status)
        call LIS_verify(status, "SMAP AC72 root zone soil moisture number of bins in the CDF: not defined")
     endif
  enddo

  do n=1,LIS_rc%nnest
     SMAP_AC72_rzmc_struc(n)%cdf_read_mon = .false.
     call ESMF_ConfigFindLabel(LIS_config, "SMAP AC72 root zone CDF read option:", rc=status)
     call ESMF_ConfigGetAttribute(LIS_config, SMAP_AC72_rzmc_struc(n)%cdf_read_opt, rc=status)
     call LIS_verify(status, "SMAP AC72 root zone CDF read option: not defined")
  enddo

  do n=1,LIS_rc%nnest
     call ESMF_AttributeSet(OBS_State(n),"Data Update Status", .false., rc=status); call LIS_verify(status)
     call ESMF_AttributeSet(OBS_State(n),"Data Update Time", -99.0, rc=status); call LIS_verify(status)
     call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status", .false., rc=status); call LIS_verify(status)
     call ESMF_AttributeSet(OBS_State(n),"Number Of Observations", LIS_rc%obs_ngrid(k), rc=status); call LIS_verify(status)
  enddo
  write(LIS_logunit,*)'[INFO] read  SMAP AC72 root zone soil moisture data specifications'

  ! Create Obs field container (1 var)
  do n=1,LIS_rc%nnest
     write(temp,'(i2.2)') 1
     read(temp,'(2a1)') vid
     obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_obsvecGrid(n,k), &
                     name="Observation"//vid(1)//vid(2), rc=status)
     call LIS_verify(status)

     ! Open obs attributes file
     write(LIS_logunit,*) '[INFO] Opening attributes for observations ', trim(LIS_rc%obsattribfile(k))
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
     read(ftn,*)
     read(ftn,*) LIS_rc%nobtypes(k)
     read(ftn,*)

     allocate(vname(LIS_rc%nobtypes(k)))
     allocate(varmax(LIS_rc%nobtypes(k)))
     allocate(varmin(LIS_rc%nobtypes(k)))
     do i=1,LIS_rc%nobtypes(k)
        read(ftn,'(a40)') vname(i)
        read(ftn,*) varmin(i), varmax(i)
        write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
     enddo
     call LIS_releaseUnitNumber(ftn)

     allocate(ssdev(LIS_rc%obs_ngrid(k)))

     if (trim(LIS_rc%perturb_obs(k)).ne."none") then
        ! Perturbation field (same logical size)
        allocate(obs_pert%vname(1)); allocate(obs_pert%perttype(1))
        allocate(obs_pert%ssdev(1)); allocate(obs_pert%stdmax(1))
        allocate(obs_pert%zeromean(1))
        allocate(obs_pert%tcorr(1)); allocate(obs_pert%xcorr(1)); allocate(obs_pert%ycorr(1))
        allocate(obs_pert%ccorr(1,1))

        call LIS_readPertAttributes(1, LIS_rc%obspertAttribfile(k), obs_pert)

        ssdev = obs_pert%ssdev(1)
        SMAP_AC72_rzmc_struc(n)%ssdev_inp = obs_pert%ssdev(1)

        pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec, grid=LIS_obsEnsOnGrid(n,k), &
                          name="Observation"//vid(1)//vid(2), rc=status)
        call LIS_verify(status)

        call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status); call LIS_verify(status)
        obs_temp(:,:) = 0.0

        call ESMF_AttributeSet(pertField(n),"Perturbation Type", obs_pert%perttype(1), rc=status); call LIS_verify(status)

        if (LIS_rc%obs_ngrid(k).gt.0) then
           call ESMF_AttributeSet(pertField(n),"Standard Deviation", ssdev, itemCount=LIS_rc%obs_ngrid(k), rc=status)
           call LIS_verify(status)
        endif

        call ESMF_AttributeSet(pertField(n),"Std Normal Max", obs_pert%stdmax(1), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean", obs_pert%zeromean(1), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale", obs_pert%tcorr(1), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(pertField(n),"X Correlation Scale", obs_pert%xcorr(1), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(pertField(n),"Y Correlation Scale", obs_pert%ycorr(1), rc=status); call LIS_verify(status)
        call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength", obs_pert%ccorr(1,:), itemCount=1, rc=status)
     endif

     deallocate(vname, varmax, varmin, ssdev)
  enddo

  write(LIS_logunit,*) '[INFO] Created States for SMAP AC72 RZMC observations'

  ! Allocate per-nest storage for projected obs and times (obs grid size)
  do n=1,LIS_rc%nnest
     SMAP_AC72_rzmc_struc(n)%nc = LIS_rc%obs_lnc(k)
     SMAP_AC72_rzmc_struc(n)%nr = LIS_rc%obs_lnr(k)
     allocate(SMAP_AC72_rzmc_struc(n)%smobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
     allocate(SMAP_AC72_rzmc_struc(n)%smtime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
     SMAP_AC72_rzmc_struc(n)%smtime = -1.0
  enddo

  ! CDF allocation (unchanged)
  do n=1,LIS_rc%nnest
     allocate(ssdev(LIS_rc%obs_ngrid(k)))
     ssdev = 0.0

     if (trim(LIS_rc%dascaloption(k)).eq."CDF matching") then
        call LIS_getCDFattributes(k, SMAP_AC72_rzmc_struc(n)%modelcdffile, SMAP_AC72_rzmc_struc(n)%ntimes, ngrid)
        if (SMAP_AC72_rzmc_struc(n)%cdf_read_opt.eq.0) then
           allocate(SMAP_AC72_rzmc_struc(n)%model_mu(LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes))
           allocate(SMAP_AC72_rzmc_struc(n)%model_sigma(LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),   SMAP_AC72_rzmc_struc(n)%ntimes))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),SMAP_AC72_rzmc_struc(n)%ntimes))
           allocate(SMAP_AC72_rzmc_struc(n)%model_xrange(LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_xrange(  LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%model_cdf(   LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_cdf(     LIS_rc%obs_ngrid(k), SMAP_AC72_rzmc_struc(n)%ntimes, SMAP_AC72_rzmc_struc(n)%nbins))
           ! read all months once
           call LIS_readMeanSigmaData(n,k, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
                 SMAP_AC72_rzmc_struc(n)%modelcdffile, "SoilMoist", SMAP_AC72_rzmc_struc(n)%model_mu, SMAP_AC72_rzmc_struc(n)%model_sigma)
           call LIS_readMeanSigmaData(n,k, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
                 SMAP_AC72_rzmc_struc(n)%obscdffile, "SoilMoist", SMAP_AC72_rzmc_struc(n)%obs_mu, SMAP_AC72_rzmc_struc(n)%obs_sigma)
           call LIS_readCDFdata(n,k, SMAP_AC72_rzmc_struc(n)%nbins, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
                 SMAP_AC72_rzmc_struc(n)%modelcdffile, "SoilMoist", SMAP_AC72_rzmc_struc(n)%model_xrange, SMAP_AC72_rzmc_struc(n)%model_cdf)
           call LIS_readCDFdata(n,k, SMAP_AC72_rzmc_struc(n)%nbins, SMAP_AC72_rzmc_struc(n)%ntimes, LIS_rc%obs_ngrid(k), &
                 SMAP_AC72_rzmc_struc(n)%obscdffile, "SoilMoist", SMAP_AC72_rzmc_struc(n)%obs_xrange, SMAP_AC72_rzmc_struc(n)%obs_cdf)
        else
           allocate(SMAP_AC72_rzmc_struc(n)%model_mu(LIS_rc%obs_ngrid(k),1))
           allocate(SMAP_AC72_rzmc_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),1))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),1))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),1))
           allocate(SMAP_AC72_rzmc_struc(n)%model_xrange(LIS_rc%obs_ngrid(k),1, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_xrange(  LIS_rc%obs_ngrid(k),1, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%model_cdf(   LIS_rc%obs_ngrid(k),1, SMAP_AC72_rzmc_struc(n)%nbins))
           allocate(SMAP_AC72_rzmc_struc(n)%obs_cdf(     LIS_rc%obs_ngrid(k),1, SMAP_AC72_rzmc_struc(n)%nbins))
           ! monthly reading done in the reader
        endif
     endif
     deallocate(ssdev)
  enddo

  ! IMPORTANT CHANGE: Do NOT build source grid here.
  ! The reader will detect the NetCDF lat/lon and call neighbor_interp_input_withgrid
  ! to produce rlat/rlon/n11 at first successful read.
  do n=1,LIS_rc%nnest
     call LIS_registerAlarm("SMAP_AC72 read alarm", 86400.0, 86400.0)
     SMAP_AC72_rzmc_struc(n)%startMode = .true.
     call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status); call LIS_verify(status)
     if (trim(LIS_rc%perturb_obs(k)).ne."none") then
        call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status); call LIS_verify(status)
     endif
  enddo
end subroutine SMAP_AC72_rzmc_setup

end module SMAP_AC72_rzmc_Mod

