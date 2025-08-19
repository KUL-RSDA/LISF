!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: CustomNcReader_Mod
!
! !DESCRIPTION:
!   This module contains a custom base reader for global netCDF products on a
!   regular latitude longitude grid stored in netCDF files, with all necessary
!   flags already applied beforehand during preprocessing.
!
!   For creating a concrete reader for a specific dataset, this reader has to
!   be subclassed in a separate file, and provide the following methods:
!   - <readername>_setup(k, OBS_State, OBS_Pert_State)
!   - <readername>_read(n, k, OBS_State, OBS_Pert_State)
!   - <readername>_write(n, k, OBS_State, OBS_Pert_State)
!   Typically, the subclassed readers should just call the respective methods
!   of the base reader. An example is given in CustomLAI_Mod.F90.
!
!   The custom base reader provides the following options for lis.config:
!
!   Custom <varname> data directory:
!       Path to directory that contains image netCDF files.
!   Custom <varname> netCDF variable name:
!       Name of the variable in the netCDF files.
!   Custom <varname> netCDF name prefix:
!       Prefix for the netCDF image file name. The image file names must follow
!       the pattern <prefix>_YYYY_MM_DD.nc
!   Custom <varname> spatial resolution
!       Spatial resolution of the product. It is assumed that the grid starts
!       at -180 + res/2 longitude and -90 + res/2 latitude.
!   Custom <varname> lat max:
!      Maximum latitude value in case data has been resampled and cropped to a
!      subdomain. (Optional)
!   Custom <varname> lat min:
!      Minimum latitude value in case data has been resampled and cropped to a
!      subdomain. (Optional)
!   Custom <varname> lon max:
!      Maximum longitude value in case data has been resampled and cropped to a
!      subdomain. (Optional)
!   Custom <varname> lon min:
!      Minimum longitude value in case data has been resampled and cropped to a
!      subdomain. (Optional)
!   Custom <varname> assimilate at local time:
!       Whether to perform the assimilation at local time (Optional, default:
!       0).
!   Custom <varname> assimilation time (hour):
!   Custom <varname> assimilation time (minute):
!       Time when to assimilate for local time DA, default is 0:00
!   Custom <varname> observation perturbation option:
!       One of "scalar", "spatial", "full".
!       - "scalar": a scalar is read from the observation perturbation
!         attributes file and used for the complete domain and period (but might
!         be rescaled according to the rescaling rules)
!       - "spatial": Observation perturbations are read from the mean
!         observation uncertainty file, observation perturbation attributes are
!         ignored. Additionally, the rescaling rules are applied.
!       - "full": Each day a new observation uncertainty is read from the image
!         files, and scaled with the product of the scalar factor in the
!         observations perturbations attributes file and the mean uncertainty
!         map. Additionally, the rescaling rules are applied.
!   Custom <varname> mean observation uncertainty file:
!       File containting uncertainty map for observation perturbation option
!       "spatial" or "full".
!   Custom <varname> mean observation uncertainty variable name:
!       Variable name of uncertainty map for observation perturbation option
!       "spatial" or "full".
!   Custom <varname> observation uncertainty variable name:
!       Name of the uncertainty variable that is read from the image files.
!   Custom <varname> use scaled standard devation model:
!       Whether to rescale the observation perturbations if observations are
!       rescaled.
!   Data assimilation scaling strategy:
!       Options are "none", "CDF matching", "Normal deviate scaling", "seasonal", "seasonal multiplicative"
!   Custom <varname> model scaling file:
!       Path to the model CDF/mean file (optional, only required if scaling is
!       applied)
!   Custom <varname> observation scaling file:
!       Path to the observation CDF/mean file (optional, only required if
!       scaling is applied)
!   Custom <varname> varname in model scaling file:
!       Prefix of the names of the variables in the model scaling files,
!       e.g. <prefix>_mu, <prefix>_sigma, ... (optional, only required if
!       scaling is applied)
!   Custom <varname> varname in observation scaling file:
!       Prefix of the names of the variables in the observation scaling files,
!       e.g. <prefix>_mu, <prefix>_sigma, ... (optional, only required if
!       scaling is applied)
!   Custom <varname> number of bins in the CDF:
!       Number of bins in the CDFs. (optional, only required if CDF scaling is
!       applied)
!
!
!
! !REVISION HISTORY:
!  02 Mar 2022    Samuel Scherrer; initial reader based on MODIS LAI reader
!
#include "LIS_misc.h"
module CustomNcReader_Mod
    ! !USES:
    use ESMF
    use map_utils
    use, intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: CustomNcReader_setup, CustomNcReader_rescale_with_seasonal_scaling,&
         CustomNcReader_updateSsdev, read_CustomNetCDF, write_CustomNetCDF
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    !EOP
    type, public:: CustomNcReader_dec

        ! This first set of attributes needs to be set in the
        ! <readername>_setup subroutine
        character*100          :: varname
        ! should be the id that is set in LIS_pluginIndices
        character*100          :: obsid
        real                   :: max_value
        real                   :: min_value
        real                   :: qcmax_value
        real                   :: qcmin_value

        integer                :: nc
        integer                :: nr
        integer                :: mi
        real                   :: gridDesci(50)
        real*8                 :: time1, time2
        integer                :: fnd
        integer                :: useSsdevScal
        logical                :: mult_scaling
        character*100          :: nc_varname
        character*100          :: nc_prefix
        real*8                 :: spatialres
        real*8                 :: dlat, dlon
        real*8                 :: latmax, latmin, lonmax, lonmin
        real,    allocatable   :: rlat(:)
        real,    allocatable   :: rlon(:)
        integer, allocatable   :: n11(:)
        integer, allocatable   :: n12(:)
        integer, allocatable   :: n21(:)
        integer, allocatable   :: n22(:)
        real,    allocatable   :: w11(:)
        real,    allocatable   :: w12(:)
        real,    allocatable   :: w21(:)
        real,    allocatable   :: w22(:)
        integer                :: lt_assim
        integer                :: da_hr, da_mn

        integer                :: obs_pert_option
        integer                :: ssdev_opt
        real                   :: ssdev_inp
        real,    allocatable   :: ssdev_inp_field(:)
        character*256          :: mean_obs_unc_file
        character*100          :: mean_obs_unc_varname
        character*100          :: obs_unc_varname
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

        real,     allocatable      :: daobs(:,:)
        real,     allocatable      :: daobs_unc(:,:)
        real,     allocatable      :: datime(:,:)

    end type CustomNcReader_dec

contains

    !BOP
    !
    ! !ROUTINE: CustomNcReader_setup
    ! \label{CustomNcReader_setup}
    !
    ! !INTERFACE:
    subroutine CustomNcReader_setup(reader_struc, k, OBS_State, OBS_Pert_State)
        ! !USES:
        use ESMF
        use LIS_coreMod
        use LIS_timeMgrMod
        use LIS_historyMod
        use LIS_dataAssimMod
        use LIS_perturbMod
        use LIS_DAobservationsMod
        use LIS_logmod

        implicit none

        ! !ARGUMENTS:
        type(CustomNcReader_dec)  :: reader_struc(LIS_rc%nnest)
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)
        !
        ! !DESCRIPTION:
        !
        !   This routine completes the runtime initializations and
        !   creation of data structures required for handling Custom NcReader data.
        !
        !   The arguments are:
        !   \begin{description}
        !    \item[k] number of observation state
        !    \item[OBS\_State]   observation state
        !    \item[OBS\_Pert\_State] observation perturbations state
        !   \end{description}
        !EOP
        integer                :: n,i
        integer                :: ftn
        integer                :: status
        type(ESMF_Field)       :: obsField(LIS_rc%nnest)
        type(ESMF_Field)       :: pertField(LIS_rc%nnest)
        type(ESMF_ArraySpec)   :: intarrspec, realarrspec
        type(ESMF_ArraySpec)   :: pertArrSpec
        character*100          :: laiobsdir
        character*100          :: temp
        real, parameter        :: minssdev =0.001
        real,  allocatable     :: ssdev(:)
        character*1            :: vid(2)
        type(pert_dec_type)    :: obs_pert
        real, pointer          :: obs_temp(:,:)
        character*40, allocatable  ::  vname(:)
        real        , allocatable  ::  varmin(:)
        real        , allocatable  ::  varmax(:)
        character*100          :: modelscalingfile(LIS_rc%nnest)
        character*100          :: obsscalingfile(LIS_rc%nnest)
        character*100          :: modelscalingvarname(LIS_rc%nnest)
        character*100          :: obsscalingvarname(LIS_rc%nnest)
        character*20           :: obs_pert_option(LIS_rc%nnest)
        integer                :: c,r,t
        integer                :: ngrid
        integer                :: timeidx
        character*100          :: varname

        ! varname is the same for all nests
        varname = reader_struc(1)%varname

        call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        !------------------------------------------------------------
        ! File options
        !------------------------------------------------------------
        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" data directory:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,laiobsdir,&
                 rc=status)
            call LIS_verify(status, &
                 "Custom "//trim(varname)//" data directory: is missing")

            call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
                 laiobsdir, rc=status)
            call LIS_verify(status)
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" netCDF variable name:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%nc_varname,&
                 rc=status)
            call LIS_verify(status,&
                 "Custom "//trim(varname)//" netCDF variable name: is missing")
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" netCDF name prefix:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%nc_prefix,&
                 rc=status)
            call LIS_verify(status, &
                 "Custom "//trim(varname)//" netCDF name prefix: is missing")
        enddo

        call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" spatial resolution:",&
             rc=status)
        do n=1,LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%spatialres,&
                 rc=status)
            call LIS_verify(status, "Custom "//trim(varname)//" spatial resolution: is missing")
        enddo

        do n=1,LIS_rc%nnest
            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" lat max:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%latmax = 90 - 0.5 * reader_struc(n)%spatialres
            else
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%latmax,&
                     rc=status)
            endif

            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" lat min:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%latmin = -90 + 0.5 * reader_struc(n)%spatialres
            else
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%latmin,&
                     rc=status)
            endif
            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" lon max:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%lonmax = 180 - 0.5 * reader_struc(n)%spatialres
            else
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%lonmax,&
                     rc=status)
            endif
            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" lon min:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%lonmin = -180 + 0.5 * reader_struc(n)%spatialres
            else
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%lonmin,&
                     rc=status)
            endif

            reader_struc(n)%dlat = reader_struc(n)%spatialres
            reader_struc(n)%dlon = reader_struc(n)%spatialres
            reader_struc(n)%nr = nint((reader_struc(n)%latmax - reader_struc(n)%latmin)&
                                       / reader_struc(n)%spatialres) + 1
            reader_struc(n)%nc = nint((reader_struc(n)%lonmax - reader_struc(n)%lonmin)&
                                       / reader_struc(n)%spatialres) + 1
        enddo

        !------------------------------------------------------------
        ! DA options
        !------------------------------------------------------------

        do n=1,LIS_rc%nnest
            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" assimilate at local time:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%lt_assim = 0
            else
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%lt_assim,&
                     rc=status)
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" assimilation time (hour):",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(reader_struc(n)%lt_assim.ne.0) then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%da_hr, &
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" assimilation time (hour): not defined")
            else
                reader_struc(n)%da_hr = 0.0
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" assimilation time (minute):",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(reader_struc(n)%lt_assim.ne.0) then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%da_mn, &
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" assimilation time (minute): not defined")
            else
                reader_struc(n)%da_mn = 0.0
            endif
        enddo

        do n=1,LIS_rc%nnest
            call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" observation perturbation option:",&
                 rc=status)
            if (status .ne. 0) then
                reader_struc(n)%obs_pert_option = 0  ! scalar
            else
                call ESMF_ConfigGetAttribute(LIS_config,obs_pert_option(n),&
                     rc=status)
                if (obs_pert_option(n).eq."scalar") then
                    reader_struc(n)%obs_pert_option = 0
                else if (obs_pert_option(n).eq."spatial") then
                    reader_struc(n)%obs_pert_option = 1
                else if (obs_pert_option(n).eq."full") then
                    reader_struc(n)%obs_pert_option = 2
                else
                    write(LIS_logunit,*)&
                         "[ERROR] Custom "//trim(varname)//" observation perturbation option"&
                         //" must be one of 'scalar', 'spatial', or 'full'"
                    call LIS_endrun
                endif
            endif
        enddo


        call ESMF_ConfigFindLabel(LIS_config,"Custom "//trim(varname)//" mean observation uncertainty file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if (reader_struc(n)%obs_pert_option .ne. 0) then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%mean_obs_unc_file,&
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" mean observation uncertainty file: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" mean observation uncertainty variable name:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(reader_struc(n)%obs_pert_option .ne. 0) then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%mean_obs_unc_varname, &
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" mean observation uncertainty variable name: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" observation uncertainty variable name:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(reader_struc(n)%obs_pert_option .eq. 2) then  ! option "full"
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%obs_unc_varname, &
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" observation uncertainty variable name: not defined")
            endif
        enddo

        !------------------------------------------------------------
        ! Options for scaling
        !------------------------------------------------------------

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" use scaled standard deviation model:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%useSsdevScal, &
                     rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" use scaled standard deviation model: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" model scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_ConfigGetAttribute(LIS_config,modelscalingfile(n),rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" model scaling file: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" observation scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_ConfigGetAttribute(LIS_config,obsscalingfile(n),rc=status)
                call LIS_verify(status,&
                     "Custom "//trim(varname)//" observation scaling file: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" varname in model scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_ConfigGetAttribute(LIS_config,modelscalingvarname(n),rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" varname in model scaling file: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config,&
             "Custom "//trim(varname)//" varname in observation scaling file:",&
             rc=status)
        do n=1,LIS_rc%nnest
            if(LIS_rc%dascaloption(k).ne."none") then
                call ESMF_ConfigGetAttribute(LIS_config,obsscalingvarname(n),rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" varname in observation scaling file: not defined")
            endif
        enddo

        call ESMF_ConfigFindLabel(LIS_config, &
             "Custom "//trim(varname)//" number of bins in the CDF:", rc=status)
        do n=1, LIS_rc%nnest
            if(LIS_rc%dascaloption(k).eq."CDF matching"&
                 .or.LIS_rc%dascaloption(k).eq."unsafe CDF matching"&
                 .or.LIS_rc%dascaloption(k).eq."Normal deviate scaling") then
                call ESMF_ConfigGetAttribute(LIS_config,reader_struc(n)%nbins, rc=status)
                call LIS_verify(status, &
                     "Custom "//trim(varname)//" number of bins in the CDF: not defined")
            endif
        enddo


        do n=1,LIS_rc%nnest
            call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
                 .false., rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
                 -99.0, rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
                 .false., rc=status)
            call LIS_verify(status)

            call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
                 LIS_rc%obs_ngrid(k),rc=status)
            call LIS_verify(status)

        enddo

        do n=1, LIS_rc%nnest
            allocate(reader_struc(n)%daobs(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
            allocate(reader_struc(n)%datime(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
            if (reader_struc(n)%obs_pert_option.eq.2) then
                allocate(reader_struc(n)%daobs_unc(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
            endif
        enddo

        write(LIS_logunit,*)&
             "[INFO] read Custom "//trim(varname)//" data specifications"

        !------------------------------------------------------------
        ! Read grid information and resample
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            if(LIS_rc%lis_obs_map_proj(k).ne."latlon") then
                write(LIS_logunit,*)&
                     '[ERROR] The Custom reader module only works with latlon projection'
                call LIS_endrun
            endif

            reader_struc(n)%gridDesci(1) = 0  ! regular lat-lon grid
            reader_struc(n)%gridDesci(2) = reader_struc(n)%nc
            reader_struc(n)%gridDesci(3) = reader_struc(n)%nr
            reader_struc(n)%gridDesci(4) = reader_struc(n)%latmax
            reader_struc(n)%gridDesci(5) = reader_struc(n)%lonmin
            reader_struc(n)%gridDesci(6) = 128
            reader_struc(n)%gridDesci(7) = reader_struc(n)%latmin
            reader_struc(n)%gridDesci(8) = reader_struc(n)%lonmax
            reader_struc(n)%gridDesci(9) = reader_struc(n)%dlat
            reader_struc(n)%gridDesci(10) = reader_struc(n)%dlon
            reader_struc(n)%gridDesci(20) = 64

            reader_struc(n)%mi = reader_struc(n)%nc*reader_struc(n)%nr

            !-----------------------------------------------------------------------------
            !   Use interpolation if LIS is running finer than native resolution.
            !-----------------------------------------------------------------------------
            if(LIS_rc%obs_gridDesc(k,10).le.reader_struc(n)%dlon) then

                allocate(reader_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
                allocate(reader_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))

                write(LIS_logunit,*)&
                     "[INFO] create interpolation input for Custom "//trim(varname)//""

                call bilinear_interp_input_withgrid(reader_struc(n)%gridDesci(:), &
                     LIS_rc%obs_gridDesc(k,:),&
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
                     reader_struc(n)%rlat, reader_struc(n)%rlon,&
                     reader_struc(n)%n11, reader_struc(n)%n12, &
                     reader_struc(n)%n21, reader_struc(n)%n22, &
                     reader_struc(n)%w11, reader_struc(n)%w12, &
                     reader_struc(n)%w21, reader_struc(n)%w22)
            else

                allocate(reader_struc(n)%n11(&
                     reader_struc(n)%nc*reader_struc(n)%nr))

                write(LIS_logunit,*)&
                     "[INFO] create upscaling input for Custom "//trim(varname)//""

                call upscaleByAveraging_input(reader_struc(n)%gridDesci(:),&
                     LIS_rc%obs_gridDesc(k,:),&
                     reader_struc(n)%nc*reader_struc(n)%nr, &
                     LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), reader_struc(n)%n11)

                write(LIS_logunit,*)&
                     "[INFO] finished creating upscaling input for Custom "//trim(varname)//""
            endif

        enddo

        !------------------------------------------------------------
        ! Read perturbation settings
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            write(unit=temp,fmt='(i2.2)') 1
            read(unit=temp,fmt='(2a1)') vid

            obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
                 grid=LIS_obsVecGrid(n,k),&
                 name="Observation"//vid(1)//vid(2), rc=status)
            call LIS_verify(status)

            !Perturbations State
            write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
                 trim(LIS_rc%obsattribfile(k))
            ftn = LIS_getNextUnitNumber()
            open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
            read(ftn,*)
            read(ftn,*) LIS_rc%nobtypes(k)
            read(ftn,*)

            allocate(vname(LIS_rc%nobtypes(k)))
            allocate(varmax(LIS_rc%nobtypes(k)))
            allocate(varmin(LIS_rc%nobtypes(k)))

            do i=1,LIS_rc%nobtypes(k)
                read(ftn,fmt='(a40)') vname(i)
                read(ftn,*) varmin(i),varmax(i)
                write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
            enddo
            call LIS_releaseUnitNumber(ftn)

            allocate(ssdev(LIS_rc%obs_ngrid(k)))

            if(trim(LIS_rc%perturb_obs(k)).ne."none") then
                allocate(obs_pert%vname(1))
                allocate(obs_pert%perttype(1))
                allocate(obs_pert%ssdev(1))
                allocate(obs_pert%stdmax(1))
                allocate(obs_pert%zeromean(1))
                allocate(obs_pert%tcorr(1))
                allocate(obs_pert%xcorr(1))
                allocate(obs_pert%ycorr(1))
                allocate(obs_pert%ccorr(1,1))

                if (reader_struc(n)%obs_pert_option.eq.0) then ! option "scalar"
                    call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
                         obs_pert)
                    reader_struc(n)%ssdev_inp = obs_pert%ssdev(1)
                    ! Set obs err to be uniform (will be rescaled later for each grid point).
                    reader_struc(n)%ssdev_inp = obs_pert%ssdev(1)
                    ssdev = obs_pert%ssdev(1)
                    write(LIS_logunit,*)&
                         "[INFO] observation perturbation size for "//trim(varname)//":",&
                         reader_struc(n)%ssdev_inp
                else if (reader_struc(n)%obs_pert_option.eq.1) then ! option "spatial"
                    call CustomNcReader_read_mean_obs_uncertainty(reader_struc, n, k, ssdev)
                    allocate(reader_struc(n)%ssdev_inp_field(LIS_rc%obs_ngrid(k)))
                    reader_struc(n)%ssdev_inp_field = ssdev
                else if (reader_struc(n)%obs_pert_option.eq.2) then ! option "full"
                    call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
                         obs_pert)
                    call CustomNcReader_read_mean_obs_uncertainty(reader_struc, n, k, ssdev)
                    allocate(reader_struc(n)%ssdev_inp_field(LIS_rc%obs_ngrid(k)))
                    ! in the "full" option, we use the spatial field for some
                    ! spatial scaling (typically the mean uncertainty over
                    ! time), and scale this by the value in the observation
                    ! perturbation attributes file, to have a tuning factor
                    do t=1,LIS_rc%obs_ngrid(k)
                        if (.not.ssdev(t).eq.LIS_rc%udef) then
                            reader_struc(n)%ssdev_inp_field(t) = obs_pert%ssdev(1) / ssdev(t)
                        else
                            reader_struc(n)%ssdev_inp_field(t) = LIS_rc%udef
                        endif
                    end do
                    ssdev = reader_struc(n)%ssdev_inp_field
                endif

                pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
                     grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
                     rc=status)
                call LIS_verify(status)

                ! initializing the perturbations to be zero
                call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
                call LIS_verify(status)
                obs_temp(:,:) = 0

                call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
                     obs_pert%perttype(1), rc=status)
                call LIS_verify(status)

                if(LIS_rc%obs_ngrid(k).gt.0) then
                    call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status)
                endif

                call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
                     obs_pert%stdmax(1), rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
                     obs_pert%zeromean(1),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
                     obs_pert%tcorr(1),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
                     obs_pert%xcorr(1),rc=status)

                call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
                     obs_pert%ycorr(1),rc=status)

                call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
                     obs_pert%ccorr(1,:),itemCount=1,rc=status)

                call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
                call LIS_verify(status)

            endif

            deallocate(vname)
            deallocate(varmax)
            deallocate(varmin)
            deallocate(ssdev)

        enddo

        !------------------------------------------------------------
        ! Initialize scaling
        !------------------------------------------------------------
        do n=1,LIS_rc%nnest

            if(LIS_rc%dascaloption(k).ne."none") then

                if(LIS_rc%dascaloption(k).eq."CDF matching"&
                     .or.LIS_rc%dascaloption(k).eq."unsafe CDF matching"&
                     .or.LIS_rc%dascaloption(k).eq."Normal deviate scaling") then
                    !------------------------------------------------------------
                    ! CDF matching
                    !------------------------------------------------------------

                    call LIS_getCDFattributes(k,modelscalingfile(n),&
                         reader_struc(n)%ntimes, ngrid)

                    allocate(reader_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))

                    allocate(reader_struc(n)%model_xrange(&
                         LIS_rc%obs_ngrid(k), reader_struc(n)%ntimes, &
                         reader_struc(n)%nbins))
                    allocate(reader_struc(n)%obs_xrange(&
                         LIS_rc%obs_ngrid(k), reader_struc(n)%ntimes, &
                         reader_struc(n)%nbins))
                    allocate(reader_struc(n)%model_cdf(&
                         LIS_rc%obs_ngrid(k), reader_struc(n)%ntimes, &
                         reader_struc(n)%nbins))
                    allocate(reader_struc(n)%obs_cdf(&
                         LIS_rc%obs_ngrid(k), reader_struc(n)%ntimes, &
                         reader_struc(n)%nbins))

                    !----------------------------------------------------------------------------
                    ! Read the model and observation CDF data
                    !----------------------------------------------------------------------------

                    ! model mean sigma
                    call LIS_readMeanSigmaData(n,k,&
                         reader_struc(n)%ntimes, &
                         LIS_rc%obs_ngrid(k), &
                         modelscalingfile(n), &
                         modelscalingvarname(n),&
                         reader_struc(n)%model_mu,&
                         reader_struc(n)%model_sigma)

                    ! observation mean sigma
                    call LIS_readMeanSigmaData(n,k,&
                         reader_struc(n)%ntimes, &
                         LIS_rc%obs_ngrid(k), &
                         obsscalingfile(n), &
                         obsscalingvarname(n),&
                         reader_struc(n)%obs_mu,&
                         reader_struc(n)%obs_sigma)

                    ! model CDF
                    call LIS_readCDFdata(n,k,&
                         reader_struc(n)%nbins,&
                         reader_struc(n)%ntimes, &
                         LIS_rc%obs_ngrid(k), &
                         modelscalingfile(n), &
                         modelscalingvarname(n),&
                         reader_struc(n)%model_xrange,&
                         reader_struc(n)%model_cdf)

                    ! observation CDF
                    call LIS_readCDFdata(n,k,&
                         reader_struc(n)%nbins,&
                         reader_struc(n)%ntimes, &
                         LIS_rc%obs_ngrid(k), &
                         obsscalingfile(n), &
                         obsscalingvarname(n),&
                         reader_struc(n)%obs_xrange,&
                         reader_struc(n)%obs_cdf)

                    !------------------------------------------------------------
                    ! end CDF matching
                    !------------------------------------------------------------

                elseif (LIS_rc%dascaloption(k).eq."seasonal"&
                     .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative") then

                    !------------------------------------------------------------
                    ! seasonal scaling
                    !------------------------------------------------------------

                    reader_struc(n)%mult_scaling = &
                         LIS_rc%dascaloption(k).eq."seasonal multiplicative"

                    reader_struc(n)%ntimes = 366

                    allocate(reader_struc(n)%model_mu(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%model_sigma(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%obs_mu(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))
                    allocate(reader_struc(n)%obs_sigma(LIS_rc%obs_ngrid(k),&
                         reader_struc(n)%ntimes))

                    call CustomNcReader_readSeasonalScalingData(n,k,&
                         reader_struc(n)%ntimes, &
                         modelscalingfile(n), &
                         modelscalingvarname(n),&
                         reader_struc(n)%model_mu,&
                         reader_struc(n)%model_sigma)

                    call CustomNcReader_readSeasonalScalingData(n,k,&
                         reader_struc(n)%ntimes, &
                         obsscalingfile(n), &
                         obsscalingvarname(n),&
                         reader_struc(n)%obs_mu,&
                         reader_struc(n)%obs_sigma)

                    if (reader_struc(n)%mult_scaling) then
                        ! When doing the multiplicative scaling, the standard
                        ! deviation in the file is not used, and the
                        ! perturbation standard deviations should be corrected
                        ! by the ratio of means instead of the ratio of sigmas.
                        reader_struc(n)%model_sigma = reader_struc(n)%model_mu
                        reader_struc(n)%obs_sigma = reader_struc(n)%obs_mu
                    endif

                endif

            endif ! dascaloption .ne. "none"
        enddo


        do n=1,LIS_rc%nnest
            call LIS_registerAlarm("Custom "//trim(varname)//" read alarm",&
                 86400.0, 86400.0)

            call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
            call LIS_verify(status, "Adding observation field failed")

        enddo

    end subroutine CustomNcReader_setup

    !BOP
    ! !ROUTINE: read_CustomNetCDF
    ! \label{read_CustomNetCDF}
    !
    ! !INTERFACE:
    subroutine read_CustomNetCDF(reader_struc, n, k, OBS_State, OBS_Pert_State)
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

        implicit none
        ! !ARGUMENTS:
        type(CustomNcReader_dec) :: reader_struc(LIS_rc%nnest)
        integer, intent(in) :: n
        integer, intent(in) :: k
        type(ESMF_State)    :: OBS_State
        type(ESMF_State)    :: OBS_Pert_State
        !
        ! !DESCRIPTION:
        !
        !  reads the Custom observations from NETCDF files.

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


        integer                :: status
        character*100          :: obsdir
        character*300          :: fname
        integer                :: cyr, cmo, cda, chr,cmn,css,cdoy
        real                   :: wt1, wt2,ts
        integer                :: count
        real*8                 :: time
        logical                :: alarmCheck, file_exists
        integer                :: c,r,i,j,p,t
        real,          pointer :: obsl(:)
        type(ESMF_Field)       :: varfield
        integer                :: gid(LIS_rc%obs_ngrid(k))
        integer                :: assimflag(LIS_rc%obs_ngrid(k))
        logical                :: data_update
        logical                :: data_upd_flag(LIS_npes)
        logical                :: data_upd_flag_local
        logical                :: data_upd
        real                   :: observations(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        real                   :: obs_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
        real                   :: obs_unscaled(LIS_rc%obs_ngrid(k))
        real, allocatable      :: observations_unc(:) !(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        real, allocatable      :: obs_unc_current(:,:) !(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
        real, allocatable      :: obs_unc_ngrid(:) !(LIS_rc%obs_ngrid(k))
        integer                :: fnd
        real                   :: timenow
        real                   :: ssdev(LIS_rc%obs_ngrid(k))
        type(ESMF_Field)       :: pertfield
        integer                :: timeidx
        real                   :: localdatime, gmtdatime, lon, dt


        call ESMF_AttributeGet(OBS_State,"Data Directory",&
             obsdir, rc=status)
        call LIS_verify(status)
        call ESMF_AttributeGet(OBS_State,"Data Update Status",&
             data_update, rc=status)
        call LIS_verify(status)

        data_upd = .false.
        if (reader_struc(n)%obs_pert_option.eq.2) then
            allocate(observations_unc(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
            allocate(obs_unc_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)))
            allocate(obs_unc_ngrid(LIS_rc%obs_ngrid(k)))
            observations_unc = LIS_rc%udef
            obs_unc_current = LIS_rc%udef
            obs_unc_ngrid = LIS_rc%udef
        endif
        obs_unscaled = LIS_rc%udef
        obs_current = LIS_rc%udef

        ! Read the data from the file at 0UTC and store it
        alarmCheck = LIS_isAlarmRinging(LIS_rc, "Custom "&
             //trim(reader_struc(n)%varname)//" read alarm")
        if(alarmCheck) then

            reader_struc(n)%daobs = LIS_rc%udef
            reader_struc(n)%datime = -1
            observations = LIS_rc%udef
            observations_unc = LIS_rc%udef

            call create_CustomNetCDF_filename(reader_struc(n)%nc_prefix,&
                 obsdir, LIS_rc%yr, LIS_rc%mo, LIS_rc%da, fname)

            inquire(file=fname,exist=file_exists)
            if(file_exists) then
                call read_CustomNetCDF_data(reader_struc, n,k, fname,observations, observations_unc)
            else
                write(LIS_logunit,*) '[WARN] Missing observation file: ',trim(fname)
            endif

            ! set daobs and datime
            reader_struc(n)%daobs  = LIS_rc%udef
            if (reader_struc(n)%obs_pert_option.eq.2) then
                reader_struc(n)%daobs_unc  = LIS_rc%udef
                do r=1,LIS_rc%obs_lnr(k)
                    do c=1,LIS_rc%obs_lnc(k)
                        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                            if(observations(c+(r-1)*LIS_rc%obs_lnc(k)).gt.0) then
                                reader_struc(n)%daobs(c,r) = &
                                     observations(c+(r-1)*LIS_rc%obs_lnc(k))
                                reader_struc(n)%daobs_unc(c, r) = &
                                     observations_unc(c+(r-1)*LIS_rc%obs_lnc(k))
                                lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                                if (reader_struc(n)%lt_assim.ne.0) then
                                    localdatime = reader_struc(n)%da_hr * 3600.0 + reader_struc(n)%da_mn * 60.0
                                    gmtdatime = localdatime - 240 * lon
                                    if (gmtdatime.lt.0) gmtdatime = gmtdatime + 86400
                                    if (gmtdatime.ge.86400) gmtdatime = gmtdatime - 86400
                                    reader_struc(n)%datime(c, r) = gmtdatime
                                else
                                    reader_struc(n)%datime(c, r) = 0
                                endif
                            endif
                        endif
                    enddo
                enddo
            else ! obs_pert_option != 2
                do r=1,LIS_rc%obs_lnr(k)
                    do c=1,LIS_rc%obs_lnc(k)
                        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                            if(observations(c+(r-1)*LIS_rc%obs_lnc(k)).gt.0) then
                                reader_struc(n)%daobs(c,r) = &
                                     observations(c+(r-1)*LIS_rc%obs_lnc(k))
                                lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                                if (reader_struc(n)%lt_assim.ne.0) then
                                    localdatime = reader_struc(n)%da_hr * 3600.0 + reader_struc(n)%da_mn * 60.0
                                    gmtdatime = localdatime - 240 * lon
                                    if (gmtdatime.lt.0) gmtdatime = gmtdatime + 86400
                                    if (gmtdatime.ge.86400) gmtdatime = gmtdatime - 86400
                                    reader_struc(n)%datime(c, r) = gmtdatime
                                else
                                    reader_struc(n)%datime(c, r) = 0
                                endif
                            endif
                        endif
                    enddo
                enddo
            endif ! obs_pert_option == 2
        endif ! alarm check

        call ESMF_StateGet(OBS_State,"Observation01",varfield,&
             rc=status)
        call LIS_verify(status, 'Error: StateGet Observation01')

        call ESMF_FieldGet(varfield,localDE=0,farrayPtr=obsl,rc=status)
        call LIS_verify(status, 'Error: FieldGet')

        fnd = 0

        ! write only the values that should be assimilated now into obs_current
        ! only difference between the two options is whether
        ! obs_unc_current should be set in the innermost loop
        if (reader_struc(n)%obs_pert_option.eq.2) then
            do r=1,LIS_rc%obs_lnr(k)
                do c=1,LIS_rc%obs_lnc(k)
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                        ! we assimilate if the assimiliation time is now or within the next
                        ! integration interval, therefore we have to subtract now
                        ! (LIC_rc%gmt) from the da time
                        dt = (reader_struc(n)%datime(c, r) - LIS_rc%gmt*3600.0)
                        if (dt.ge.0.and.dt.lt.LIS_rc%ts) then
                            ! if local time assimilation is off, dt == 0 at gmt == 0
                            obs_current(c, r) = reader_struc(n)%daobs(c, r)
                            obs_unc_current(c, r) = reader_struc(n)%daobs_unc(c, r)
                            if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                                obs_unscaled(LIS_obs_domain(n,k)%gindex(c,r)) = &
                                     obs_current(c,r)
                            endif
                            if(obs_current(c,r).ne.LIS_rc%udef) then
                                fnd = 1
                            endif
                        endif
                    endif
                enddo
            enddo
        else !obs_pert_option.ne.2
            do r=1,LIS_rc%obs_lnr(k)
                do c=1,LIS_rc%obs_lnc(k)
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                        ! we assimilate if the assimiliation time is now or within the next
                        ! integration interval, therefore we have to subtract
                        ! "now" (LIS_rc%gmt) from the DA time
                        dt = (reader_struc(n)%datime(c, r) - LIS_rc%gmt*3600.0)
                        if (dt.ge.0.and.dt.lt.LIS_rc%ts) then
                            ! if local time assimilation is off, dt == 0 at gmt == 0
                            obs_current(c, r) = reader_struc(n)%daobs(c, r)
                            if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                                obs_unscaled(LIS_obs_domain(n,k)%gindex(c,r)) = &
                                     obs_current(c,r)
                            endif
                            if(obs_current(c,r).ne.LIS_rc%udef) then
                                fnd = 1
                            endif
                        endif
                    endif
                enddo
            enddo
        endif !obs_pert_option.eq.2


        !-------------------------------------------------------------------------
        !  Transform data to the LSM climatology using a CDF-scaling approach
        !-------------------------------------------------------------------------

        if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then

            write(LIS_logunit,*) '[INFO] perform CDF matching'
            call LIS_rescale_with_CDF_matching(     &
                 n,k,                               &
                 reader_struc(n)%nbins,         &
                 reader_struc(n)%ntimes,        &
                 reader_struc(n)%max_value, &
                 reader_struc(n)%min_value, &
                 reader_struc(n)%model_xrange,  &
                 reader_struc(n)%obs_xrange,    &
                 reader_struc(n)%model_cdf,     &
                 reader_struc(n)%obs_cdf,       &
                 obs_current)

        elseif(LIS_rc%dascaloption(k).eq."unsafe CDF matching".and.fnd.ne.0) then

            write(LIS_logunit,*) '[INFO] perform CDF matching'
            call CustomNcReader_rescale_with_unsafe_CDF_matching(     &
                 n,k,                               &
                 reader_struc(n)%nbins,         &
                 reader_struc(n)%ntimes,        &
                 reader_struc(n)%max_value, &
                 reader_struc(n)%min_value, &
                 reader_struc(n)%model_xrange,  &
                 reader_struc(n)%obs_xrange,    &
                 reader_struc(n)%model_cdf,     &
                 reader_struc(n)%obs_cdf,       &
                 obs_current)

        elseif (LIS_rc%dascaloption(k).eq."Normal deviate scaling".and.fnd.ne.0) then

            write(LIS_logunit,*) '[INFO] perform normal deviate scaling'
            call LIS_rescale_with_normal_deviate_scaling(     &
                 n,k,                               &
                 reader_struc(n)%ntimes,        &
                 reader_struc(n)%model_mu,  &
                 reader_struc(n)%model_sigma,    &
                 reader_struc(n)%obs_mu,     &
                 reader_struc(n)%obs_sigma,       &
                 obs_current)

        elseif ((LIS_rc%dascaloption(k).eq."seasonal"&
                 .or.LIS_rc%dascaloption(k).eq."seasonal multiplicative")&
             .and.fnd.ne.0) then

            write(LIS_logunit,*) '[INFO] perform seasonal rescaling'
            call CustomNcReader_rescale_with_seasonal_scaling(&
                 n,k,&
                 nint(LIS_get_curr_calday(LIS_rc, 0)), &
                 reader_struc(n)%ntimes,        &
                 reader_struc(n)%max_value, &
                 reader_struc(n)%min_value, &
                 reader_struc(n)%mult_scaling, &
                 reader_struc(n)%model_mu,  &
                 reader_struc(n)%model_sigma,  &
                 reader_struc(n)%obs_mu,  &
                 reader_struc(n)%obs_sigma,  &
                 obs_current)

        endif

        !-------------------------------------------------------------------------
        !  End transforming
        !-------------------------------------------------------------------------

        ! write transformed data to ESMF pointer
        if (reader_struc(n)%obs_pert_option.eq.2) then
            obsl = LIS_rc%udef
            obs_unc_ngrid = LIS_rc%udef
            do r=1, LIS_rc%obs_lnr(k)
                do c=1, LIS_rc%obs_lnc(k)
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                        obsl(LIS_obs_domain(n,k)%gindex(c,r)) = obs_current(c,r)
                        obs_unc_ngrid(LIS_obs_domain(n,k)%gindex(c,r)) = obs_unc_current(c,r)
                    endif
                enddo
            enddo
        else
            obsl = LIS_rc%udef
            do r=1, LIS_rc%obs_lnr(k)
                do c=1, LIS_rc%obs_lnc(k)
                    if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                        obsl(LIS_obs_domain(n,k)%gindex(c,r)) = obs_current(c,r)
                    endif
                enddo
            enddo
        endif ! obs_pert_option.eq.2

        !-------------------------------------------------------------------------
        !  Apply LSM-based QC and screening of observations
        !-------------------------------------------------------------------------
        call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
             //trim(reader_struc(n)%obsid)//char(0),n, k,OBS_state)
        call LIS_checkForValidObs(n,k,obsl,fnd,obs_current)


        data_upd_flag_local = (fnd.ne.0)
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
                call ESMF_AttributeSet(varfield,"Grid Number",&
                     gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(varfield,"Assimilation Flag",&
                     assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                call LIS_verify(status)

                call ESMF_AttributeSet(varfield, "Unscaled Obs",&
                     obs_unscaled, itemCount=LIS_rc%obs_ngrid(k), rc=status)
                call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

            endif

            ! rescale perturbation standard deviations if rescaling of the
            ! data is performed, or if we use daily different perturbations
            ! (option 2)
            if ((LIS_rc%dascaloption(k).ne."none".and.reader_struc(n)%useSsdevScal.eq.1)&
                 .or.reader_struc(n)%obs_pert_option.eq.2) then
                call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                     rc=status)
                call LIS_verify(status, 'Error: StateGet Observation01')

                if (reader_struc(n)%obs_pert_option.eq.0) then
                    ssdev = reader_struc(n)%ssdev_inp
                else if (reader_struc(n)%obs_pert_option.eq.1) then
                    ssdev = reader_struc(n)%ssdev_inp_field
                else ! option 2, i.e. "full"
                    do t=1,LIS_rc%obs_ngrid(k)
                        if (.not.obs_unc_ngrid(t).eq.LIS_rc%udef&
                            .and. .not. reader_struc(n)%ssdev_inp_field(t).eq.LIS_rc%udef) then
                            ssdev(t) = reader_struc(n)%ssdev_inp_field(t) * obs_unc_ngrid(t)
                        else
                            ssdev(t) = LIS_rc%udef
                        endif
                    end do
                endif

                if (LIS_rc%dascaloption(k).ne."none".and.reader_struc(n)%useSsdevScal.eq.1) then
                    timeidx = CustomNcReader_timeidx(reader_struc(n)%ntimes)
                    call CustomNcReader_updateSsdev(k,&
                         reader_struc(n)%obs_sigma(:, timeidx),&
                         reader_struc(n)%model_sigma(:, timeidx),&
                         ssdev)
                endif

                if(LIS_rc%obs_ngrid(k).gt.0) then
                    call ESMF_AttributeSet(pertfield,"Standard Deviation",&
                         ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
                    call LIS_verify(status)
                endif
            endif

        else ! no data update
            call ESMF_AttributeSet(OBS_State,"Data Update Status",&
                 .false., rc=status)
            call LIS_verify(status)
        endif
        if (reader_struc(n)%obs_pert_option.eq.2) then
            deallocate(observations_unc)
            deallocate(obs_unc_current)
            deallocate(obs_unc_ngrid)
        endif

    end subroutine read_CustomNetCDF

    !BOP
    !
    ! !ROUTINE: read_CustomNetCDF_data
    ! \label{read_CustomNetCDF_data}
    !
    ! !INTERFACE:
    subroutine read_CustomNetCDF_data(reader_struc, n, k, fname, obs_ip, obs_unc_ip)
        !
        ! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_coreMod,  only : LIS_rc, LIS_domain
        use LIS_logMod
        use LIS_timeMgrMod

        implicit none
        !
        ! !INPUT PARAMETERS:
        !
        type(CustomNcReader_dec)      :: reader_struc(LIS_rc%nnest)
        integer                       :: n
        integer                       :: k
        character (len=*)             :: fname
        real, intent(inout)           :: obs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        real, intent(inout)           :: obs_unc_ip(:) !(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


        ! !OUTPUT PARAMETERS:
        !
        !
        ! !DESCRIPTION:
        !  This subroutine reads the Custom netCDF file and applies the data
        !  quality flags to filter the data.
        !
        !  The arguments are:
        !  \begin{description}
        !  \item[n]            index of the nest
        !  \item[k]            number of observation state
        !  \item[k]            number of observation state
        !  \item[fname]        name of the netCDF file
        !  \item[observations\_ip]   Observation data processed to the LIS domain
        !  \end{description}
        !
        !
        !EOP

        integer                 :: lat_offset, lon_offset
        real                    :: observation(reader_struc(n)%nc,reader_struc(n)%nr)
        real                    :: obs_in(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1               :: obs_b_in(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1               :: obs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        real, allocatable       :: observation_unc(:,:) !(reader_struc(n)%nc,reader_struc(n)%nr)
        real, allocatable       :: obs_unc_in(:) !(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1, allocatable  :: obs_unc_b_in(:) !(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1, allocatable  :: obs_unc_b_ip(:) !(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        integer                 :: c,r,t
        integer                 :: nid, lid, latsize_file, lonsize_file
        integer                 :: lat(reader_struc(n)%nr)
        integer                 :: obsid, flagid
        integer                 :: ios
        integer                 :: numvalidobs

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] read_CustomNetCDF requires NETCDF"
        call LIS_endrun
#else

        if (reader_struc(n)%obs_pert_option.eq.2) then
            allocate(observation_unc(reader_struc(n)%nc,reader_struc(n)%nr))
            allocate(obs_unc_in(reader_struc(n)%nc*reader_struc(n)%nr))
            allocate(obs_unc_b_in(reader_struc(n)%nc*reader_struc(n)%nr))
            allocate(obs_unc_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
            obs_unc_b_in = .false.
        endif
        obs_b_in = .false.

        lat_offset = 1  ! no offset
        lon_offset = 1

        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        call LIS_verify(nf90_inq_dimid(nid, "lat",lid), &
             "Error nf90_inq_dimid: lat")
        call LIS_verify(nf90_inquire_dimension(nid, lid, len=latsize_file), &
             "Error nf90_inquire_dimension:lat")

        call LIS_verify(nf90_inq_dimid(nid, "lon",lid), &
             "Error nf90_inq_dimid: lon")
        call LIS_verify(nf90_inquire_dimension(nid, lid, len=lonsize_file), &
             "Error nf90_inquire_dimension:lon")

        ! check if the grid has the right size
        if (latsize_file.ne.reader_struc(n)%nr.or.lonsize_file.ne.reader_struc(n)%nc) then
            write(LIS_logunit,*) "[ERR] Dimension sizes not consistent with expectations"
            write(LIS_logunit,*) 'lat size on file: ',latsize_file
            write(LIS_logunit,*) 'lat size expected: ',reader_struc(n)%nr
            write(LIS_logunit,*) 'lon size on file: ',lonsize_file
            write(LIS_logunit,*) 'lon size expected: ',reader_struc(n)%nc
            call LIS_endrun
        endif

        ! make sure that latitude axis is inverted
        ios = nf90_inq_varid(nid, "lat", lid)
        call LIS_verify(ios, 'Error nf90_inq_varid: lat')
        ios = nf90_get_var(nid, lid, lat)
        call LIS_verify(ios, 'Error nf90_inquire_variable: lat')
        if (lat(1) < lat(latsize_file)) then
            write(LIS_logunit,*) "[ERR] Reader expects inverted latitude coordinate"
            call LIS_endrun
        endif


        ! read main variable
        ios = nf90_inq_varid(nid, trim(reader_struc(n)%nc_varname), obsid)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//reader_struc(n)%nc_varname)
        ios = nf90_get_var(nid, obsid, observation, &
             start=(/lon_offset,lat_offset/), &
             count=(/reader_struc(n)%nc,reader_struc(n)%nr/))
        call LIS_verify(ios, 'Error nf90_get_var: '//reader_struc(n)%nc_varname)

        if (reader_struc(n)%obs_pert_option.eq.2) then
            ! read uncertainty variable
            ios = nf90_inq_varid(nid, trim(reader_struc(n)%obs_unc_varname), obsid)
            call LIS_verify(ios, 'Error nf90_inq_varid: '//reader_struc(n)%obs_unc_varname)
            ios = nf90_get_var(nid, obsid, observation_unc, &
                 start=(/lon_offset,lat_offset/), &
                 count=(/reader_struc(n)%nc,reader_struc(n)%nr/))
            call LIS_verify(ios, 'Error nf90_get_var: '//reader_struc(n)%obs_unc_varname)
        endif

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))


        ! the data is already read into 'observation', but we have to replace
        ! NaNs/invalid values with LIS_rc%udef
        if (reader_struc(n)%obs_pert_option.eq.2) then
            do r=1, reader_struc(n)%nr
                do c=1, reader_struc(n)%nc
                    if (isnan(observation(c, r)) &
                         .or.isnan(observation_unc(c, r)) &
                         .or.(observation(c, r) < reader_struc(n)%qcmin_value) &
                         .or.(observation(c, r) > reader_struc(n)%qcmax_value) &
                         .or.(observation_unc(c, r) <= 0)) then
                        observation(c, r) = LIS_rc%udef
                        observation_unc(c, r) = LIS_rc%udef
                    endif
                    ! fill obs_in and obs_b_in, which are required further on
                    obs_in(c+(r-1)*reader_struc(n)%nc) = observation(c,r)
                    obs_b_in(c+(r-1)*reader_struc(n)%nc) = observation(c, r).ne.LIS_rc%udef
                    obs_unc_in(c+(r-1)*reader_struc(n)%nc) = observation_unc(c,r)
                    obs_unc_b_in(c+(r-1)*reader_struc(n)%nc) = observation(c, r).ne.LIS_rc%udef
                end do
            end do
        else !obs_pert_option.ne.2
            do r=1, reader_struc(n)%nr
                do c=1, reader_struc(n)%nc
                    if (isnan(observation(c, r)) &
                         .or.(observation(c, r) < reader_struc(n)%qcmin_value) &
                         .or.(observation(c, r) > reader_struc(n)%qcmax_value)) then
                        observation(c, r) = LIS_rc%udef
                    endif
                    ! fill obs_in and obs_b_in, which are required further on
                    obs_in(c+(r-1)*reader_struc(n)%nc) = observation(c,r)
                    obs_b_in(c+(r-1)*reader_struc(n)%nc) = observation(c, r).ne.LIS_rc%udef
                end do
            end do

        endif !obs_pert_option.eq.2

        call CustomNcReader_interp_data(reader_struc, n, k, obs_in, obs_b_in, obs_ip, obs_b_ip, fname)
        if (reader_struc(n)%obs_pert_option.eq.2) then
            call CustomNcReader_interp_data(reader_struc, n, k, obs_unc_in, obs_unc_b_in, obs_unc_ip, obs_unc_b_ip, fname)
        endif

        if (reader_struc(n)%obs_pert_option.eq.2) then
            deallocate(observation_unc)
            deallocate(obs_unc_in)
            deallocate(obs_unc_b_in)
            deallocate(obs_unc_b_ip)
        endif


        write(LIS_logunit,*) '[INFO] Finished reading ',trim(fname)
#endif
    end subroutine read_CustomNetCDF_data


    subroutine CustomNcReader_interp_data(reader_struc, n, k, data_in, data_b_in, data_ip, data_b_ip, fname)
        use LIS_coreMod,  only : LIS_rc, LIS_domain
        use LIS_logMod

        type(CustomNcReader_dec)   :: reader_struc(LIS_rc%nnest)
        integer,   intent(in)         :: n
        integer,   intent(in)         :: k
        real, intent(inout)        :: data_in(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1, intent(inout)   :: data_b_in(reader_struc(n)%nc*reader_struc(n)%nr)
        real, intent(inout)        :: data_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        logical*1, intent(inout)   :: data_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        character(len=*), intent(in) :: fname

        integer                 :: ios

        if(LIS_rc%obs_gridDesc(k,10).le.reader_struc(n)%dlon) then
            write(LIS_logunit,*) '[INFO] interpolating Custom ',&
                 trim(reader_struc(n)%varname), ' ',&
                 trim(fname)
            !--------------------------------------------------------------------------
            ! Interpolate to the LIS running domain if model has finer resolution
            ! than observations
            !--------------------------------------------------------------------------
            call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
                 data_b_in, data_in, data_b_ip, data_ip, &
                 reader_struc(n)%nc*reader_struc(n)%nr, &
                 LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
                 reader_struc(n)%rlat,reader_struc(n)%rlon,&
                 reader_struc(n)%w11,reader_struc(n)%w12,&
                 reader_struc(n)%w21,reader_struc(n)%w22,&
                 reader_struc(n)%n11,reader_struc(n)%n12,&
                 reader_struc(n)%n21,reader_struc(n)%n22,LIS_rc%udef,ios)
         else
            write(LIS_logunit,*) '[INFO] upscaling Custom',&
                 trim(reader_struc(n)%varname),&
                 trim(fname)
            !--------------------------------------------------------------------------
            ! Upscale to the LIS running domain if model has coarser resolution
            ! than observations
            !--------------------------------------------------------------------------
            call upscaleByAveraging(reader_struc(n)%nc*reader_struc(n)%nr,&
                 LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
                 LIS_rc%udef, reader_struc(n)%n11,&
                 data_b_in,data_in, data_b_ip, data_ip)
        endif
    end subroutine CustomNcReader_interp_data



    !BOP
    ! !ROUTINE: create_CustomNetCDF_filename
    ! \label{create_CustomNetCDF_filename}
    !
    ! !INTERFACE:
    subroutine create_CustomNetCDF_filename(prefix, ndir, year, month, day, filename)
        ! !USES:

        implicit none
        ! !ARGUMENTS:
        character (len=*), intent(in)     :: prefix
        character (len=*), intent(in)     :: ndir
        integer, intent(in)               :: year, month, day
        character(len=*), intent(out)     :: filename
        !
        ! !DESCRIPTION:
        !  This subroutine creates the netCDF filename
        !  based on the time and date
        !
        !  The arguments are:
        !  \begin{description}
        !  \item[varname] variable name of the netCDF variable
        !  \item[ndir] name of the data directory
        !  \item[year]  current year
        !  \item[month]  current month
        !  \item[day]  current day
        !  \item[filename] Generated filename
        !  \end{description}
        !
        !EOP

        !  The naming scheme is
        !    <varname>_<year>_<month>_<day>.nc


        character(len=4) :: yearstr
        character(len=2) :: monthstr
        character(len=2) :: daystr

        write(unit=yearstr, fmt='(i4.4)') year
        write(unit=monthstr, fmt='(i2.2)') month
        write(unit=daystr, fmt='(i2.2)') day


        filename = trim(ndir)//'/'//&
             trim(prefix)//'_'//yearstr//'_'//monthstr//'_'//daystr//'.nc'

    end subroutine create_CustomNetCDF_filename

    ! !ROUTINE: write_CustomNetCDF
    ! \label{write_CustomNetCDF}
    !
    ! !INTERFACE:
    subroutine write_CustomNetCDF(n, k, OBS_State)
    ! !USES:
      use ESMF
      use LIS_coreMod
      use LIS_logMod
      use LIS_fileIOMod
      use LIS_historyMod
      use LIS_DAobservationsMod

      implicit none

    ! !ARGUMENTS:

      integer,     intent(in)  :: n
      integer,     intent(in)  :: k
      type(ESMF_State)         :: OBS_State
    !
    ! !DESCRIPTION:
    !
    ! writes the transformed (interpolated/upscaled/reprojected)
    ! observations to a file
    !
    !EOP
      type(ESMF_Field)         :: obsField
      logical                  :: data_update
      real, pointer            :: observations(:)
      character*100            :: obsname
      integer                  :: ftn
      integer                  :: status

      call ESMF_AttributeGet(OBS_State, "Data Update Status", &
           data_update, rc=status)
      call LIS_verify(status)

      if(data_update) then

         call ESMF_StateGet(OBS_State, "Observation01",obsField, &
              rc=status)
         call LIS_verify(status)

         call ESMF_FieldGet(obsField, localDE=0, farrayPtr=observations, rc=status)
         call LIS_verify(status)

         if(LIS_masterproc) then
            ftn = LIS_getNextUnitNumber()
            call Custom_obsname(n,k,obsname)

            call LIS_create_output_directory('DAOBS')
            open(ftn,file=trim(obsname), form='unformatted')
         endif

         call LIS_writevar_gridded_obs(ftn,n,k,observations)

         if(LIS_masterproc) then
            call LIS_releaseUnitNumber(ftn)
         endif

      endif

    end subroutine write_CustomNetCDF

    !BOP
    ! !ROUTINE: Custom_obsname
    ! \label{Custom_obsname}
    !
    ! !INTERFACE:
    subroutine Custom_obsname(n,k,obsname)
    ! !USES:
      use LIS_coreMod, only : LIS_rc

    ! !ARGUMENTS:
      integer               :: n
      integer               :: k
      character(len=*)      :: obsname
    !
    ! !DESCRIPTION:
    !
    !  writes the assimilated observation to a file.

    !
    !  The arguments are:
    !  \begin{description}
    !  \item[n] index of the nest
    !  \item[k] number of observation state
    !  \item[obsname] name of the observation
    !  \end{description}
    !

    !EOP

      character(len=12) :: cdate1
      character(len=12) :: cdate
      character(len=10) :: cda

      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LIS_rc%yr, LIS_rc%mo, &
           LIS_rc%da, LIS_rc%hr,LIS_rc%mn

      write(unit=cda, fmt='(a2,i2.2)') '.a',k
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n

      obsname = trim(LIS_rc%odir)//'/DAOBS/'//cdate1(1:6)//&
           '/LISDAOBS_'//cdate1// &
           trim(cda)//trim(cdate)//'.1gs4r'

    end subroutine Custom_obsname


    !BOP
    !
    ! !ROUTINE: CustomNcReader_rescale_with_seasonal_scaling
    ! \label{CustomNcReader_rescale_with_seasonal_scaling}
    !
    ! !INTERFACE:
    subroutine CustomNcReader_rescale_with_seasonal_scaling(&
         n,             &
         k,             &
         timeidx,        &
         ntimes,         &
         max_obs_value, &
         min_obs_value, &
         multiplicative, &
         model_mu,    &
         model_sigma,    &
         obs_mu,       &
         obs_sigma,       &
         obs_value)

        use LIS_coreMod, only: LIS_rc
        use LIS_DAobservationsMod, only: LIS_obs_domain

        implicit none
        !
        ! !ARGUMENTS:
        integer, intent(in)      :: n
        integer, intent(in)      :: k
        integer, intent(in)      :: ntimes
        integer, intent(in)      :: timeidx
        real, intent(in)         :: max_obs_value
        real, intent(in)         :: min_obs_value
        logical, intent(in)      :: multiplicative
        real, intent(in)         :: model_mu(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: model_sigma(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: obs_mu(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(in)         :: obs_sigma(LIS_rc%obs_ngrid(k),ntimes)
        real, intent(inout)      :: obs_value(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
        !
        ! !DESCRIPTION:
        !
        !   This routine rescales the input observation data by subtracting the
        !   observation mean for the given time index and adding the model mean
        !   for the same time index

        !  The arguments are:
        !  \begin{description}
        !  \item[n]               index of the nest
        !  \item[k]               index of observation state
        !  \item[timeidx]         time index for mean value to use, .e.g. month
        !                         or day
        !  \item[ntimes]          number of times in the model and obs mean files
        !  \item[max\_obs\_value] maximum allowable value of observation
        !  \item[min\_obs\_value] minimum allowable value of observation
        !  \item[model\_mu]       mean values of the model OL run
        !  \item[model\_sigma]    std.dev. values of the model OL run
        !  \item[obs\_mu]         mean values of the observations
        !  \item[obs\_sigma]      std.dev. values of the observations
        !  \item[obs\_value]      observation value to be rescaled.
        ! \end{description}
        !EOP

        integer             :: grididx,i
        integer             :: binval
        integer             :: col,row
        real                :: cdf_obsval
        real                :: obs_tmp

        do grididx=1,LIS_rc%obs_ngrid(k)

            col = LIS_obs_domain(n,k)%col(grididx)
            row = LIS_obs_domain(n,k)%row(grididx)

            if((.not. isnan(obs_value(col,row))).and.obs_value(col,row).ne.-9999.0 &
                 .and. obs_mu(grididx, timeidx).ne.-9999.0 &
                 .and. obs_sigma(grididx, timeidx).gt.epsilon(0.0)) then

                if (multiplicative.and.obs_mu(grididx, timeidx).ne.epsilon(0.0)) then
                    obs_tmp = obs_value(col, row) / obs_mu(grididx, timeidx)&
                         * model_mu(grididx, timeidx)
                else if ((.not.multiplicative)&
                        .and. obs_sigma(grididx, timeidx).gt.(epsilon(0.0))) then
                    ! transform to z-score
                    obs_tmp = (obs_value(col,row) - obs_mu(grididx,timeidx)) &
                         / obs_sigma(grididx,timeidx)
                    ! transform z-score to model space
                    obs_tmp = obs_tmp * model_sigma(grididx,timeidx)&
                         + model_mu(grididx, timeidx)
                else
                    obs_tmp = LIS_rc%udef
                endif

                if (obs_tmp < min_obs_value .or. obs_tmp > max_obs_value) then
                    obs_tmp = LIS_rc%udef
                endif

                obs_value(col, row) = obs_tmp
            else
                obs_value(col,row) = LIS_rc%udef
            endif
        enddo
    end subroutine CustomNcReader_rescale_with_seasonal_scaling

    !BOP
    ! !ROUTINE: CustomNcReader_readSeasonalScalingData
    ! \label{CustomNcReader_readSeasonalScalingData}
    !
    ! !INTERFACE:
    subroutine CustomNcReader_readSeasonalScalingData(&
         n, k, ntimes, filename, varname, mu, sigma)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_coreMod, only: LIS_rc
        use LIS_logMod, only: LIS_logunit, LIS_verify, LIS_endrun
        use LIS_DAobservationsMod, only: LIS_convertObsVarToLocalSpace

        implicit none
        ! !ARGUMENTS:
        integer,   intent(in)         :: n
        integer,   intent(in)         :: k
        integer,   intent(in)         :: ntimes
        character(len=*), intent(in)  :: filename
        character(len=*), intent(in)  :: varname
        real, intent(inout)           :: mu(:,:)
        real, intent(inout)           :: sigma(:,:)
        !
        ! !DESCRIPTION:
        !  This routine reads the input seasonal mean file.
        !
        !  The arguments are:
        !  \begin{description}
        !  \item[n]             index of the nest
        !  \item[k]             index of observation state
        !  \item[ngrid]         length of ngrid dimension
        !  \item[ntimes]        length of ntimes dimension
        !  \item[filename]      name of the CDF file
        !  \item[varname]       name of the variable being extracted.
        !  \item[mu]            mean values
        !  \item[sigma]         std.dev values
        ! \end{description}
        !EOP
        integer                  :: j
        integer                  :: nlevsId, gId
        integer                  :: ngrid_file, nlevs_file
        integer                  :: muid, sigmaid
        real, allocatable        :: mu_file(:,:,:), sigma_file(:,:,:)
        integer                  :: nid

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] read_CustomNetCDF requires NETCDF"
        call LIS_endrun
#else
        write(LIS_logunit,*) "[INFO] Reading mean from seasonal seasonal file ",trim(filename)
        call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=nid),"failed to open file "//trim(filename))

        call LIS_verify(nf90_inq_dimid(nid,trim(varname)//"_levels",nlevsId), &
             "nf90_inq_dimid failed for "//trim(varname)//"_levels")

        call LIS_verify(nf90_inquire_dimension(nid,nlevsId, len=nlevs_file),&
             "nf90_inquire_dimension failed for nlevsId")

        call LIS_verify(nf90_inq_dimid(nid, "ngrid",gId), &
             "Error nf90_inq_dimid: ngrid")

        call LIS_verify(nf90_inquire_dimension(nid, gId, len=ngrid_file), &
             "Error nf90_inquire_dimension:ngrid")

        if (ngrid_file /= LIS_rc%obs_glbngrid_red(k)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(filename)//" not consistent "//&
                 "with expected ngrid: ", ngrid_file,&
                 " instead of ",LIS_rc%obs_glbngrid_red(k)
            call LIS_endrun
        endif

        ! dimension order is flipped compared to netCDF, because in Fortran the
        ! first dimension changes fastest
        allocate(mu_file(ngrid_file,ntimes,nlevs_file))
        allocate(sigma_file(ngrid_file,ntimes,nlevs_file))

        call LIS_verify(nf90_inq_varid(nid,trim(varname)//"_mu",muid),&
             "nf90_inq_varid failed for for "//trim(varname)//"_mu")

        call LIS_verify(nf90_get_var(nid,muid,mu_file),&
             "nf90_get_var failed for "//trim(varname)//"_mu")

        call LIS_verify(nf90_inq_varid(nid,trim(varname)//"_sigma",sigmaid),&
             "nf90_inq_varid failed for for "//trim(varname)//"_sigma")

        call LIS_verify(nf90_get_var(nid,sigmaid,sigma_file),&
             "nf90_get_var failed for "//trim(varname)//"_sigma")


        if(LIS_rc%obs_ngrid(k).gt.0) then
            do j=1,ntimes
                call LIS_convertObsVarToLocalSpace(n,k,mu_file(:,j,1), mu(:,j))
                call LIS_convertObsVarToLocalSpace(n,k,sigma_file(:,j,1), sigma(:,j))
            enddo
        endif

        deallocate(mu_file)
        deallocate(sigma_file)

        call LIS_verify(nf90_close(nid),&
             "failed to close file "//trim(filename))
        write(LIS_logunit,*)&
             "[INFO] Successfully read mean and seasonal scaling file ",&
             trim(filename)
#endif
    end subroutine CustomNcReader_readSeasonalScalingData

    !BOP
    ! !ROUTINE: CustomNcReader_readSsdevData
    ! \label{CustomNcReader_readSsdevData}
    !
    ! !INTERFACE:
    subroutine CustomNcReader_readSsdevData(&
         n, k, filename, varname, ssdev)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_coreMod, only: LIS_rc
        use LIS_logMod, only: LIS_logunit, LIS_verify, LIS_endrun
        use LIS_DAobservationsMod, only: LIS_convertObsVarToLocalSpace

        implicit none
        ! !ARGUMENTS:
        integer,   intent(in)         :: n
        integer,   intent(in)         :: k
        character(len=*), intent(in)  :: filename
        character(len=*), intent(in)  :: varname
        real, intent(inout)           :: ssdev(:)
        !
        ! !DESCRIPTION:
        !  This routine reads the mean uncertainty file
        !
        !  The arguments are:
        !  \begin{description}
        !  \item[n]             index of the nest
        !  \item[k]             index of observation state
        !  \item[ngrid]         length of ngrid dimension
        !  \item[filename]      name of the CDF file
        !  \item[varname]       name of the variable being extracted.
        !  \item[ssdev]         observation perturbation values
        ! \end{description}
        !EOP
        integer                  :: j
        integer                  :: gId
        integer                  :: ngrid_file
        integer                  :: ssdevid
        real, allocatable        :: ssdev_file(:)
        integer                  :: nid

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] read_CustomNetCDF requires NETCDF"
        call LIS_endrun
#else
        write(LIS_logunit,*) "[INFO] Reading ssdev from file ",trim(filename)
        call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
             ncid=nid),"failed to open file "//trim(filename))

        call LIS_verify(nf90_inq_dimid(nid, "ngrid",gId), &
             "Error nf90_inq_dimid: ngrid")

        call LIS_verify(nf90_inquire_dimension(nid, gId, len=ngrid_file), &
             "Error nf90_inquire_dimension:ngrid")

        if (ngrid_file /= LIS_rc%obs_glbngrid_red(k)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(filename)//" not consistent "//&
                 "with expected ngrid: ", ngrid_file,&
                 " instead of ",LIS_rc%obs_glbngrid_red(k)
            call LIS_endrun
        endif

        allocate(ssdev_file(ngrid_file))

        call LIS_verify(nf90_inq_varid(nid,trim(varname),ssdevid),&
             "nf90_inq_varid failed for for "//trim(varname))

        call LIS_verify(nf90_get_var(nid,ssdevid,ssdev_file),&
             "nf90_get_var failed for "//trim(varname))

        if(LIS_rc%obs_ngrid(k).gt.0) then
            call LIS_convertObsVarToLocalSpace(n,k,ssdev_file(:), ssdev(:))
        endif

        deallocate(ssdev_file)

        call LIS_verify(nf90_close(nid),&
             "failed to close file "//trim(filename))
        write(LIS_logunit,*)&
             "[INFO] Successfully read ssdev file ",&
             trim(filename)
#endif
    end subroutine CustomNcReader_readSsdevData

    subroutine CustomNcReader_updateSsdev(k, obs_sigma, model_sigma, ssdev)
        use LIS_coreMod
        implicit none

        integer, intent(in)     :: k
        real, intent(in)        :: obs_sigma(:), model_sigma(:)
        real, intent(inout)     :: ssdev(:)

        integer                 :: grididx
        real, parameter         ::  minssdev =0.001


        do grididx=1,LIS_rc%obs_ngrid(k)
            if(obs_sigma(grididx).ne.LIS_rc%udef) then
                if(obs_sigma(grididx).ne.0) then
                    ssdev(grididx) = ssdev(grididx)&
                         * model_sigma(grididx)&
                         / obs_sigma(grididx)
                endif

                if(ssdev(grididx).lt.minssdev) then
                    ssdev(grididx) = minssdev
                endif
            endif
        enddo

    end subroutine CustomNcReader_updateSsdev

    function CustomNcReader_timeidx(ntimes)
        use LIS_coreMod,  only : LIS_rc
        use LIS_logMod
        use LIS_timeMgrMod
        implicit none

        integer, intent(in) :: ntimes
        integer :: CustomNcReader_timeidx

        if (ntimes == 1) then
            CustomNcReader_timeidx = 1
        else if (ntimes == 12) then
            CustomNcReader_timeidx = LIS_rc%mo
        else if (ntimes == 366) then
            CustomNcReader_timeidx = nint(LIS_get_curr_calday(LIS_rc, 0))
        else
            write(LIS_logunit,*)&
                 "[ERR] Unexpected number of times in scaling file:", ntimes
            call LIS_endrun
        endif

    end function CustomNcReader_timeidx

    !BOP
    !
    ! !ROUTINE: read_CustomNetCDF_data
    ! \label{read_CustomNetCDF_data}
    !
    ! !INTERFACE:
    subroutine CustomNcReader_read_mean_obs_uncertainty(reader_struc, n, k, unc_ngrid)
        !
        ! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_coreMod,  only : LIS_rc, LIS_domain
        use LIS_logMod
        use LIS_timeMgrMod
        use LIS_DAobservationsMod, only: LIS_obs_domain

        implicit none
        !
        ! !INPUT PARAMETERS:
        !
        type(CustomNcReader_dec)      :: reader_struc(LIS_rc%nnest)
        integer,   intent(in)         :: n
        integer,   intent(in)         :: k
        real, intent(inout)           :: unc_ngrid(:)

        integer                 :: lat_offset, lon_offset
        real                    :: uncertainty(reader_struc(n)%nc,reader_struc(n)%nr)
        real                    :: unc_in(reader_struc(n)%nc*reader_struc(n)%nr)
        real                    :: unc_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        logical*1               :: unc_b_in(reader_struc(n)%nc*reader_struc(n)%nr)
        logical*1               :: unc_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
        integer                 :: c,r,t
        integer                 :: nid, lid, latsize_file, lonsize_file
        integer                 :: lat(reader_struc(n)%nr)
        integer                 :: obsid, flagid
        integer                 :: ios
        integer                 :: numvalidobs
        character*256           :: fname
        character*100           :: varname

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] read_CustomNetCDF requires NETCDF"
        call LIS_endrun
#else

        fname = reader_struc(n)%mean_obs_unc_file
        varname = reader_struc(n)%mean_obs_unc_varname

        unc_b_in = .false.

        lat_offset = 1  ! no offset
        lon_offset = 1

        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '//trim(fname))

        call LIS_verify(nf90_inq_dimid(nid, "lat",lid), &
             "Error nf90_inq_dimid: lat")
        call LIS_verify(nf90_inquire_dimension(nid, lid, len=latsize_file), &
             "Error nf90_inquire_dimension:lat")

        call LIS_verify(nf90_inq_dimid(nid, "lon",lid), &
             "Error nf90_inq_dimid: lon")
        call LIS_verify(nf90_inquire_dimension(nid, lid, len=lonsize_file), &
             "Error nf90_inquire_dimension:lon")

        ! check if the grid has the right size
        if (latsize_file.ne.reader_struc(n)%nr.or.lonsize_file.ne.reader_struc(n)%nc) then
            write(LIS_logunit,*) "[ERR] Dimension sizes not consistent with expectations"
            write(LIS_logunit,*) 'lat size on file: ',latsize_file
            write(LIS_logunit,*) 'lat size expected: ',reader_struc(n)%nr
            write(LIS_logunit,*) 'lon size on file: ',lonsize_file
            write(LIS_logunit,*) 'lon size expected: ',reader_struc(n)%nc
            call LIS_endrun
        endif

        ! make sure that latitude axis is inverted
        ios = nf90_inq_varid(nid, "lat", lid)
        call LIS_verify(ios, 'Error nf90_inq_varid: lat')
        ios = nf90_get_var(nid, lid, lat)
        call LIS_verify(ios, 'Error nf90_inquire_variable: lat')
        if (lat(1) < lat(latsize_file)) then
            write(LIS_logunit,*) "[ERR] Reader expects inverted latitude coordinate"
            call LIS_endrun
        endif

        ! read main variable
        ios = nf90_inq_varid(nid, trim(varname), obsid)
        call LIS_verify(ios, 'Error nf90_inq_varid: '//varname)
        ios = nf90_get_var(nid, obsid, uncertainty, &
             start=(/lon_offset,lat_offset/), &
             count=(/reader_struc(n)%nc,reader_struc(n)%nr/))
        call LIS_verify(ios, 'Error nf90_get_var: '//varname)

        ios = nf90_close(ncid=nid)
        call LIS_verify(ios,'Error closing file '//trim(fname))


        ! the data is already read into 'uncertainty', but we have to replace
        ! NaNs/invalid values with LIS_rc%udef
        do r=1, reader_struc(n)%nr
            do c=1, reader_struc(n)%nc
                if (isnan(uncertainty(c, r))) then
                    uncertainty(c, r) = LIS_rc%udef
                endif
                ! fill obs_in and obs_b_in, which are required further on
                unc_in(c+(r-1)*reader_struc(n)%nc) = uncertainty(c,r)
                unc_b_in(c+(r-1)*reader_struc(n)%nc) = uncertainty(c, r).ne.LIS_rc%udef
            end do
        end do

        call CustomNcReader_interp_data(reader_struc, n, k, unc_in, unc_b_in, unc_ip, unc_b_ip, fname)
        write(LIS_logunit,*) '[INFO] Finished reading ',trim(fname)

        ! now we have to move from the (obs_lnc(k)*obs_lnr(k),) shape to
        ! (obs_ngrid(k),)
        unc_ngrid = LIS_rc%udef
        do r=1, LIS_rc%obs_lnr(k)
            do c=1, LIS_rc%obs_lnc(k)
                if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then
                    unc_ngrid(LIS_obs_domain(n,k)%gindex(c,r)) = unc_ip(c+(r-1)*LIS_rc%obs_lnc(k))
                endif
             end do
         end do
#endif
    end subroutine CustomNcReader_read_mean_obs_uncertainty

!BOP
!
! !ROUTINE: CustomNcReader_rescale_with_unsafe_CDF_matching
! \label{CustomNcReader_rescale_with_unsafe_CDF_matching}
!
! !INTERFACE:
  subroutine CustomNcReader_rescale_with_unsafe_CDF_matching(&
       n,             &
       k,             &
       nbins,         &
       ntimes,        &
       max_obs_value, &
       min_obs_value, &
       model_xrange,  &
       obs_xrange,    &
       model_cdf,     &
       obs_cdf,       &
       obs_value)

    use LIS_coreMod,  only : LIS_rc
    use LIS_DAobservationsMod, only: LIS_obs_domain
    use LIS_logmod

    implicit none
!
! !ARGUMENTS:
    integer             :: n
    integer             :: k
    integer             :: nbins
    integer             :: ntimes
    real                :: max_obs_value
    real                :: min_obs_value
    real                :: model_xrange(LIS_rc%obs_ngrid(k),ntimes, nbins)
    real                :: obs_xrange(LIS_rc%obs_ngrid(k),ntimes, nbins)
    real                :: model_cdf(LIS_rc%obs_ngrid(k),ntimes, nbins)
    real                :: obs_cdf(LIS_rc%obs_ngrid(k),ntimes, nbins)
    real                :: obs_value(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
!
! !DESCRIPTION:
!
!   This is a copy of LIS_rescale_with_CDF_matching but without the check for
!   CDF steepness.
!
!  The arguments are:
!  \begin{description}
!  \item[n]               index of the nest
!  \item[nbins]           number of bins used to compute the model and obs CDFs
!  \item[max\_obs\_value] maximum allowable value of observation
!  \item[min\_obs\_value] minimum allowable value of observation
!  \item[model\_xrange]   x-axis values corresponding to the model CDF
!  \item[obs\_xrange]     x-axis values corresponding to the obs CDF
!  \item[model\_cdf]      y-axis (CDF) values corresponding to the model CDF
!  \item[obs\_cdf]        y-axis (CDF) values corresponding to the obs CDF
!  \item[obs\_value]      observation value to be rescaled.
! \end{description}
!EOP

    real                :: model_delta(LIS_rc%obs_ngrid(k))
    real                :: obs_delta(LIS_rc%obs_ngrid(k))

    integer             :: t,i,kk
    integer             :: binval
    integer             :: col,row
    real                :: cdf_obsval
    real                :: obs_tmp, obs_in
    integer, dimension (1)      :: index_25 , index_75
    real                :: Lb_xrange, Ub_xrange, iqr_obs, iqr_model

    if(ntimes.gt.1) then
       kk = LIS_rc%mo
    else
       kk = 1
    endif
    do t=1,LIS_rc%obs_ngrid(k)
       model_delta(t) = model_xrange(t,kk,2)-model_xrange(t,kk,1)
       obs_delta(t)   = obs_xrange(t,kk,2)-obs_xrange(t,kk,1)
    enddo
    do t=1,LIS_rc%obs_ngrid(k)

        col = LIS_obs_domain(n,k)%col(t)
        row = LIS_obs_domain(n,k)%row(t)
! MN : Replaced the dynamic range with the interquartile range

!SVK dynamic range of obs is too small.
!        if((obs_xrange(t,kk,nbins)-obs_xrange(t,kk,1)).lt.&
!           (max_obs_value-min_obs_value)*0.25) then
!           obs_delta(t) = 0
!        endif

        index_25 = minloc(abs(obs_cdf(t,kk,:) - 0.25))
        index_75 = minloc(abs(obs_cdf(t,kk,:) - 0.75))
        Lb_xrange = obs_xrange(t,kk,index_25(1))
        Ub_xrange = obs_xrange(t,kk,index_75(1))
        iqr_obs = Ub_xrange - Lb_xrange

        index_25 = minloc(abs(model_cdf(t,kk,:) - 0.25))
        index_75 = minloc(abs(model_cdf(t,kk,:) - 0.75))
        Lb_xrange = model_xrange(t,kk,index_25(1))
        Ub_xrange = model_xrange(t,kk,index_75(1))
        iqr_model = Ub_xrange - Lb_xrange

        if(obs_value(col,row).ne.-9999.0) then
           obs_in = obs_value(col,row)
           if(obs_delta(t).gt.0) then
              binval = nint((obs_value(col,row)-obs_xrange(t,kk,1))/&
                   obs_delta(t))+1
              if(binval.gt.nbins) binval = nbins
              if(binval.le.0) binval = 1
              cdf_obsval = obs_cdf(t,kk,binval)
              if(cdf_obsval.gt.1.0) cdf_obsval = 1.0
              i=1
              do while((model_cdf(t,kk,i).lt.cdf_obsval).and.&
                   (i.le.nbins))
                 i = i+1
                 if(i.gt.nbins) exit
              enddo
              if(i.gt.nbins) i = i-1
              obs_tmp = model_xrange(t,kk,i)

              if(obs_tmp.gt.max_obs_value) then
!                 obs_tmp = max_obs_value
                 obs_tmp = LIS_rc%udef
              endif

              if(obs_tmp.le.min_obs_value) then
!                 obs_tmp = obs_value(col,row)
                 obs_tmp = LIS_rc%udef
              endif
              obs_value(col,row) = obs_tmp
           else
              obs_value(col,row) = LIS_rc%udef
           endif
           if(obs_value(col,row).le.min_obs_value.and.&
                obs_value(col,row).ne.-9999.0) then
              write(LIS_logunit,*) '[ERR] Problem in CDF scaling of observations in the DA instance ',k
              write(LIS_logunit,*) '[ERR] ',col,row,obs_value(col,row), obs_in
              call LIS_endrun()
           endif
        else
           obs_value(col,row) = LIS_rc%udef
        endif
     enddo
   end subroutine CustomNcReader_rescale_with_unsafe_CDF_matching

end module CustomNcReader_Mod
