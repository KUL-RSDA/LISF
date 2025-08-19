!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module VODOO_Mod
    !BOP
    !
    ! !MODULE: VODOO_Mod
    !
    ! !DESCRIPTION:
    !    This modules implements an observation operator for VOD.
    !
    !    Although it's not a real RTM it programmatically works similarly and is
    !    therefore implemented analogously to a RTM.
    !
    !    The module provides the following options for lis.config:
    !
    !    VODOO parameter file:
    !        Path to the parameter file
    !    VODOO model type:
    !        Type of the model to use, either "window" or "anom"
    !
    !
    ! !HISTORY:
    ! 02 Mar 2022: Samuel Scherrer; initial contribution based on WCMRTM
    ! 25 Aug 2022: Samuel Scherrer; support for multiple models
    ! 17 Feb 2023: Samuel Scherrer; complete rewrite with new model structure
    !
    ! !USES:

#if (defined RTMS)

    use ESMF
    use LIS_coreMod
    use LIS_RTMMod
    use LIS_logMod

    implicit none

    PRIVATE

    !-----------------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    !-----------------------------------------------------------------------------
    public :: VODOO_initialize
    public :: VODOO_f2t
    public :: VODOO_run
    public :: VODOO_output
    public :: VODOO_geometry
    !-----------------------------------------------------------------------------
    ! !PUBLIC TYPES:
    !-----------------------------------------------------------------------------
    public :: vodoo_struc
    !EOP

    ! The actual implementation of the model equations is done in the
    ! subclasses below, this is just to provide a common interface.
    ! Using instances of this type will raise an error.
    type, public ::  vodoo_type_dec
        character*256 :: parameter_fname
        !-------output------------!
        real, allocatable :: VOD(:)
        integer :: ncoef
        contains
        procedure, pass(self) :: initialize => VODOO_initialize_default
        procedure, pass(self) :: run => VODOO_run_default
        procedure, pass(self) :: isvalid => VODOO_isvalid_default
    end type vodoo_type_dec

    type, extends(vodoo_type_dec) :: vodoo_window_model
        ! coefs has dimensions (ngrid, doy, predictor)
        real, allocatable :: coefs(:, :, :)
        real, allocatable :: intercept(:, :)
        contains
            procedure, pass(self) :: initialize => VODOO_initialize_window_model
            procedure, pass(self) :: run => VODOO_run_window_model
            procedure, pass(self) :: isvalid => VODOO_isvalid_window_model
    end type vodoo_window_model

    type, extends(vodoo_type_dec) :: vodoo_anom_model
        real, allocatable :: vod_clim(:, :)
        real, allocatable :: lai_clim(:, :)
        real, allocatable :: rzsm_clim(:, :)
        real, allocatable :: coefs(:, :)
        real, allocatable :: intercept(:)
        contains
            procedure, pass(self) :: initialize => VODOO_initialize_anom_model
            procedure, pass(self) :: run => VODOO_run_anom_model
            procedure, pass(self) :: isvalid => VODOO_isvalid_anom_model
    end type vodoo_anom_model

    class(vodoo_type_dec), allocatable :: vodoo_struc(:)

    SAVE

contains
    !BOP
    !
    ! !ROUTINE: VODOO_initialize
    ! \label{VODOO_initialize}
    !
    ! !INTERFACE:
    subroutine VODOO_initialize()
        ! !DESCRIPTION:
        !
        !  This routine creates the datatypes and allocates memory for noahMP3.6-specific
        !  variables. It also invokes the routine to read the runtime specific options
        !  for noahMP3.6 from the configuration file.
        !
        !  The routines invoked are:
        !  \begin{description}
        !   \item[readVODOOcrd](\ref{readVODOOcrd}) \newline
        !
        !EOP
        implicit none

        integer :: rc, ios
        integer :: n, nid, ngrid, ngridId
        character*100 :: modeltype(LIS_rc%nnest)
        character*256 :: parameter_fname(LIS_rc%nnest)

        write(LIS_logunit,*) "[INFO] Starting VODOO setup"

        call ESMF_ConfigFindLabel(LIS_config, "VODOO model type:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, modeltype(n), rc=rc)
            call LIS_verify(rc, "VODOO model type: not defined")
        enddo
        if (modeltype(1) .eq. "window") then
            allocate(vodoo_window_model :: vodoo_struc(LIS_rc%nnest))
        elseif (modeltype(1) .eq. "anom") then
            allocate(vodoo_anom_model :: vodoo_struc(LIS_rc%nnest))
        else
            write(LIS_logunit, *)&
                 "[ERR] VODOO model type must be 'window' or 'anom'"
            call LIS_endrun
        endif

        ! read config from file
        call ESMF_ConfigFindLabel(LIS_config, "VODOO parameter file:",rc = rc)
        do n=1, LIS_rc%nnest
            call ESMF_ConfigGetAttribute(LIS_config, parameter_fname(n), rc=rc)
            call LIS_verify(rc, "VODOO parameter file: not defined")
        enddo


        do n=1,LIS_rc%nnest
            vodoo_struc(n)%parameter_fname = parameter_fname(n)
            allocate(vodoo_struc(n)%VOD(LIS_rc%npatch(n,LIS_rc%lsm_index)))

            call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")
            call add_sfc_fields(n,LIS_sfcState(n), "Root Zone Soil Moisture")
            call add_sfc_fields(n,LIS_forwardState(n),"VODOO_VOD")
        enddo


        ! read parameters/initialize model
        do n=1, LIS_rc%nnest
            call vodoo_struc(n)%initialize(n)
        enddo

        write(LIS_logunit,*) '[INFO] Finished VODOO setup'
    end subroutine VODOO_initialize

    subroutine VODOO_initialize_default(self, n)
        class(vodoo_type_dec), intent(inout) :: self
        integer, intent(in) :: n
        write(LIS_logunit,*) "[ERR] VODOO should use window or anom model"
        call LIS_endrun
    end subroutine VODOO_initialize_default


    subroutine VODOO_initialize_window_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodoo_window_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid

        self%ncoef = 2

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodoo_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! now that we have ntimes, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%coefs(self%ncoef, 366, npatch))
        allocate(self%intercept(366, npatch))

        call read_3d_coef_from_file(n, nid, self%ncoef, 366, ngrid, "coefs", &
            self%coefs, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, 366, ngrid, "intercept", &
            self%intercept, ngrid_first=.false.)
#endif
    end subroutine VODOO_initialize_window_model


    subroutine VODOO_initialize_anom_model(self, n)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        class(vodoo_anom_model), intent(inout) :: self
        integer, intent(in) :: n

        integer :: ios, nid, npatch
        integer :: ngridId, ngrid

        self%ncoef = 5

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        ! try opening the parameter file
        write(LIS_logunit,*) '[INFO] Reading ',&
            trim(self%parameter_fname)
        ios = nf90_open(path=trim(self%parameter_fname),&
            mode=NF90_NOWRITE,ncid=nid)
        call LIS_verify(ios,'Error opening file '&
            //trim(vodoo_struc(n)%parameter_fname))

        ! check if ngrid is as expected
        ios = nf90_inq_dimid(nid, "ngrid", ngridId)
        call LIS_verify(ios, "Error nf90_inq_varid: ngrid")
        ios = nf90_inquire_dimension(nid, ngridId, len=ngrid)
        call LIS_verify(ios, "Error nf90_inquire_dimension: ngrid")
        if (ngrid /= LIS_rc%glbngrid_red(n)) then
            write(LIS_logunit, *) "[ERR] ngrid in "//trim(self%parameter_fname)&
                 //" not consistent with expected ngrid: ", ngrid,&
                 " instead of ",LIS_rc%glbngrid_red(n)
            call LIS_endrun
        endif

        ! now that we have ntimes, we can allocate the coefficient arrays
        ! for the nest
        npatch = LIS_rc%npatch(n, LIS_rc%lsm_index)
        allocate(self%vod_clim(366, npatch))
        allocate(self%lai_clim(366, npatch))
        allocate(self%rzsm_clim(366, npatch))
        allocate(self%coefs(5, npatch))
        allocate(self%intercept(npatch))

        call read_2d_coef_from_file(n, nid, 366, ngrid, "vod_climatology", &
            self%vod_clim, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, 366, ngrid, "lai_climatology", &
            self%lai_clim, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, 366, ngrid, "rzsm_climatology", &
            self%rzsm_clim, ngrid_first=.false.)
        call read_2d_coef_from_file(n, nid, self%ncoef, ngrid, "coefs", &
            self%coefs, ngrid_first=.false.)
        call read_1d_coef_from_file(n, nid, ngrid, "intercept", self%intercept)
#endif
    end subroutine VODOO_initialize_anom_model


    subroutine read_1d_coef_from_file(n, nid, ngrid, varname, coef)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, ngrid
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:)

        real, allocatable :: coef_file(:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else

        allocate(coef_file(ngrid))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        call LIS_convertVarToLocalSpace(n, coef_file, coef_grid)
        call gridvar_to_patchvar(&
             n, LIS_rc%lsm_index, coef_grid, coef)
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_1d_coef_from_file

    subroutine read_2d_coef_from_file(n, nid, n1, n2, varname, coef, ngrid_first)
        ! Parameters
        ! ----------
        ! n: nest index
        ! nid: netcdf file id
        ! n1: size of dim 1
        ! n2: size of dim 2
        ! varname: name of the variable
        ! coef: array to fill (inout)
        ! ngrid_first: indicates whether ngrid is the first or last dimension
        ! on the file
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))

        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j=1,n2
                call LIS_convertVarToLocalSpace(n, coef_file(:,j), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(:, j))
            enddo
        else ! ngrid is the second dimension, we have to loop over n1
            do j=1,n1
                call LIS_convertVarToLocalSpace(n, coef_file(j, :), coef_grid)
                call gridvar_to_patchvar(&
                     n, LIS_rc%lsm_index, coef_grid, coef(j, :))
            enddo
        endif

        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_2d_coef_from_file

    subroutine read_3d_coef_from_file(n, nid, n1, n2, n3, varname, coef, &
            ngrid_first)
        ! Parameters
        ! ----------
        ! n: nest index
        ! nid: netcdf file id
        ! n1: size of dim 1
        ! n2: size of dim 2
        ! n3: size of dim 3
        ! varname: name of the variable
        ! coef: array to fill (inout)
        ! ngrid_first: indicates whether ngrid is the first or last dimension
        ! on the file
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        use netcdf
#endif
        use LIS_historyMod, only: LIS_convertVarToLocalSpace
        integer, intent(in) :: n, nid, n1, n2, n3
        character(len=*), intent(in) :: varname
        real, intent(inout) :: coef(:,:,:)
        logical, intent(in), optional :: ngrid_first

        logical :: ngrid_is_first
        real, allocatable :: coef_file(:,:,:)
        real, allocatable :: coef_grid(:)
        integer :: ios, varid, j2, j3

#if !(defined USE_NETCDF3 || defined USE_NETCDF4)
        write(LIS_logunit,*) "[ERR] VODOO requires NETCDF"
        call LIS_endrun
#else
        if (present(ngrid_first)) then
            ngrid_is_first = ngrid_first
        else
            ngrid_is_first = .true.
        endif

        allocate(coef_file(n1, n2, n3))
        ios = nf90_inq_varid(nid, trim(varname), varid)
        call LIS_verify(ios, "Error nf90_inq_varid: "//trim(varname))
        ios = nf90_get_var(nid, varid, coef_file)
        call LIS_verify(ios, "Error nf90_get_var: "//trim(varname))
        allocate(coef_grid(LIS_rc%ngrid(n)))
        if (ngrid_is_first) then
            do j2=1,n2
                do j3=1,n3
                    call LIS_convertVarToLocalSpace(n, coef_file(:,j2,j3), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(:, j2,j3))
                 enddo
            enddo
        else  ! ngrid is last
            do j2=1,n1
                do j3=1,n2
                    call LIS_convertVarToLocalSpace(n, coef_file(j2,j3, :), coef_grid)
                    call gridvar_to_patchvar(&
                         n, LIS_rc%lsm_index, coef_grid, coef(j2,j3, :))
                 enddo
            enddo
        endif
        deallocate(coef_grid)
        deallocate(coef_file)
#endif
    end subroutine read_3d_coef_from_file

    subroutine gridvar_to_patchvar(n,m,gvar,tvar)
        ! Converts a variable in local gridspace (length LIS_rc%ngrid(n))
        ! to local patch space (length LIS_rc%npatch(n,m))
        ! patch space = ensembles * ngrid

        implicit none

        integer, intent(in) :: n 
        integer, intent(in) :: m
        real, intent(in)    :: gvar(LIS_rc%ngrid(n))
        real, intent(inout) :: tvar(LIS_rc%npatch(n,m))
        integer             :: t,r,c

        do t=1,LIS_rc%npatch(n,m)
            r = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            c = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            if (LIS_domain(n)%gindex(c,r).ge.0) then
                tvar(t) = gvar(LIS_domain(n)%gindex(c,r))
            else
                tvar(t) = LIS_rc%udef
            endif
        enddo

    end subroutine gridvar_to_patchvar
    !!--------------------------------------------------------------------------------



    subroutine add_sfc_fields(n, sfcState,varname)

        implicit none

        integer            :: n
        type(ESMF_State)   :: sfcState
        character(len=*)   :: varname

        type(ESMF_Field)     :: varField
        type(ESMF_ArraySpec) :: arrspec
        integer              :: status
        real :: sum
        call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
             rc=status)
        call LIS_verify(status)

        varField = ESMF_FieldCreate(arrayspec=arrSpec, &
             grid=LIS_vecTile(n), name=trim(varname), &
             rc=status)
        call LIS_verify(status, 'Error in field_create of '//trim(varname))

        call ESMF_StateAdd(sfcState, (/varField/), rc=status)
        call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

    end subroutine add_sfc_fields

    subroutine VODOO_f2t(n)

        implicit none

        integer, intent(in)    :: n

    end subroutine VODOO_f2t

    subroutine VODOO_geometry(n)
        implicit none
        integer, intent(in)    :: n

    end subroutine VODOO_geometry

    subroutine VODOO_run(n)
        use LIS_histDataMod
        ! !USES:
        implicit none

        integer, intent(in) :: n

        integer             :: t
        integer             :: status
        integer             :: col,row
        real, pointer       :: lai(:), rzsm(:)
        real, pointer       :: vodval(:)
        real                :: smanom
        real                :: pred(vodoo_struc(n)%ncoef)

        !   map surface properties to SFC
        call getsfcvar(LIS_sfcState(n), "Leaf Area Index", lai)
        call getsfcvar(LIS_sfcState(n), "Root Zone Soil Moisture", rzsm)

        !---------------------------------------------
        ! Patch loop
        !--------------------------------------------
        do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)

            pred(1) = lai(t)
            pred(2) = rzsm(t)

            if (vodoo_struc(n)%isvalid(n, t)) then
                call vodoo_struc(n)%run(n, t, pred)
            else
                vodoo_struc(n)%VOD(t) = LIS_rc%udef
            endif

            if (vodoo_struc(n)%VOD(t).ne.LIS_rc%udef.and.vodoo_struc(n)%VOD(t).lt.-10) then
                write(LIS_logunit, *) "[WARN] VOD lower than -10"
                vodoo_struc(n)%VOD(t) = LIS_rc%udef
            endif

            call LIS_diagnoseRTMOutputVar(n, t, LIS_MOC_RTM_VOD,&
                 value=vodoo_struc(n)%VOD(t),&
                 vlevel=1,&
                 unit="-",&
                 direction="-")
        enddo

        call getsfcvar(LIS_forwardState(n),"VODOO_VOD", vodval)
        vodval = vodoo_struc(n)%VOD
    end subroutine VODOO_run

    subroutine VODOO_run_default(self, n, t, pred)
        class(vodoo_type_dec), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%ncoef)
        write(LIS_logunit,*) "[ERR] VODOO should use window or anom model"
        call LIS_endrun
    end subroutine VODOO_run_default

    subroutine VODOO_run_window_model(self, n, t, pred)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodoo_window_model), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%ncoef)

        real                :: coefs(self%ncoef)
        integer             :: doy 

        doy = VODOO_current_doy()
        coefs = self%coefs(:, doy, t)
        self%VOD(t) = sum(coefs * pred) + self%intercept(doy, t)
    end subroutine VODOO_run_window_model


    subroutine VODOO_run_anom_model(self, n, t, pred)
        use LIS_histDataMod
        ! !USES:
        implicit none

        class(vodoo_anom_model), intent(inout) :: self
        integer, intent(in) :: n, t
        real, intent(in) :: pred(self%ncoef)

        real                :: coefs(self%ncoef)
        real                :: actual_preds(self%ncoef)
        real                :: vod_clim, lai_clim, rzsm_clim
        real                :: lai, rzsm, lai_anom, rzsm_anom
        integer             :: doy 

        doy = VODOO_current_doy()

        coefs = self%coefs(:, t)
        vod_clim = self%vod_clim(doy, t)
        lai_clim = self%lai_clim(doy, t)
        rzsm_clim = self%rzsm_clim(doy, t)

        lai = pred(1)
        rzsm = pred(2)
        lai_anom = lai - lai_clim
        rzsm_anom = rzsm - rzsm_clim
        
        actual_preds(1) = lai_anom
        actual_preds(2) = rzsm_anom
        actual_preds(3) = lai * rzsm_anom
        actual_preds(4) = rzsm * rzsm_anom
        actual_preds(5) = lai * rzsm * rzsm_anom

        self%VOD(t) = sum(coefs * actual_preds) + self%intercept(t) + vod_clim
    end subroutine VODOO_run_anom_model

    function VODOO_current_doy()
        use LIS_coreMod, only : LIS_rc

        integer   :: VODOO_current_doy
        integer   :: k, days(13), ldays(13)

        data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
        data ldays /31,29,31,30,31,30,31,31,30,31,30,31,30/

        VODOO_current_doy = 0
        if((mod(LIS_rc%yr,4).eq.0.and.mod(LIS_rc%yr,100).ne.0) &     !correct for leap year
             .or.(mod(LIS_rc%yr,400).eq.0))then             !correct for y2k
            do k=1,(LIS_rc%mo-1)
                VODOO_current_doy=VODOO_current_doy+ldays(k)
            enddo
        else
            do k=1,(LIS_rc%mo-1)
                VODOO_current_doy=VODOO_current_doy+days(k)
            enddo
        endif
        VODOO_current_doy = VODOO_current_doy + LIS_rc%da
    end function VODOO_current_doy

    function VODOO_isvalid_default(self, n, t)
        class(vodoo_type_dec), intent(inout) :: self
        integer, intent(in) :: n, t
        logical             :: VODOO_isvalid_default
        VODOO_isvalid_default = .false.
        write(LIS_logunit,*) "[ERR] VODOO should use window or anom model"
        call LIS_endrun
    end function VODOO_isvalid_default

    function VODOO_isvalid_window_model(self, n, t)
        class(vodoo_window_model), intent(inout) :: self
        integer, intent(in) :: n, t
        logical             :: VODOO_isvalid_window_model

        integer             :: doy 
        real                :: coefs(self%ncoef)

        doy = VODOO_current_doy()
        coefs = self%coefs(:, doy, t)

        VODOO_isvalid_window_model = (.not.isnan(coefs(1)).and.coefs(1).ne.LIS_rc%udef)
    end function VODOO_isvalid_window_model

    function VODOO_isvalid_anom_model(self, n, t)
        class(vodoo_anom_model), intent(inout) :: self
        integer, intent(in) :: n, t
        logical             :: VODOO_isvalid_anom_model

        integer             :: doy 
        real                :: vod_clim, lai_clim, rzsm_clim
        real                :: coefs(self%ncoef)

        doy = VODOO_current_doy()
        coefs = self%coefs(:, t)
        vod_clim = self%vod_clim(doy, t)
        lai_clim = self%lai_clim(doy, t)
        rzsm_clim = self%rzsm_clim(doy, t)

        ! if any coef is NaN or undefined, all are, but climatologies might differ
        VODOO_isvalid_anom_model = (.not.isnan(coefs(1)).and.coefs(1).ne.LIS_rc%udef&
                .and..not.isnan(vod_clim).and.vod_clim.ne.LIS_rc%udef&
                .and..not.isnan(lai_clim).and.lai_clim.ne.LIS_rc%udef&
                .and..not.isnan(rzsm_clim).and.rzsm_clim.ne.LIS_rc%udef)
    end function VODOO_isvalid_anom_model


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine getsfcvar(sfcState, varname, var)
        ! !USES:

        implicit none

        type(ESMF_State)      :: sfcState
        type(ESMF_Field)      :: varField
        character(len=*)      :: varname
        real, pointer         :: var(:)
        integer               :: status

        call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
        call LIS_verify(status, "Error in StateGet: VODOO_getsfcvar "//trim(varname))
        call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
        call LIS_verify(status, "Error in FieldGet: VODOO_getsfcvar "//trim(varname))

    end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: VODOO_output
!!!! \label{VODOO_output}
!!!!
!!!! !INTERFACE:
    subroutine VODOO_output(n)
        integer, intent(in) :: n
    end subroutine VODOO_output
#endif
end module VODOO_Mod
