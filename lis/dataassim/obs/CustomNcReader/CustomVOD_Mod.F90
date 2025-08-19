module CustomVOD_Mod
    use CustomNcReader_Mod, only: CustomNcReader_dec

    implicit none

    public :: CustomVOD_setup, read_CustomVOD, write_CustomVOD
    public :: CustomVOD_struc

    ! declare public reader array
    type(CustomNcReader_dec), allocatable :: CustomVOD_struc(:)

contains

    subroutine CustomVOD_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomVOD_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomVOD_struc(n)%obsid = "Custom VOD"
            CustomVOD_struc(n)%varname = "VOD"
            CustomVOD_struc(n)%min_value = 0.0001
            CustomVOD_struc(n)%max_value = 4.0
            CustomVOD_struc(n)%qcmin_value = 0.0
            CustomVOD_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(CustomVOD_struc, k, OBS_State, OBS_Pert_State)

    end subroutine CustomVOD_setup

    subroutine read_CustomVOD(n, k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: read_CustomNetCDF
        integer, intent(in)       :: n, k
        type(ESMF_State)          :: OBS_State
        type(ESMF_State)          :: OBS_Pert_State
        call read_CustomNetCDF(CustomVOD_struc, n, k, OBS_State, OBS_Pert_State)
    end subroutine read_CustomVOD

    subroutine write_CustomVOD(n, k, OBS_State)
        use ESMF
        use CustomNcReader_Mod, only: write_CustomNetCDF
        integer,     intent(in)  :: n, k
        type(ESMF_State)         :: OBS_State
        call write_CustomNetCDF(n, k, OBS_State)
    end subroutine write_CustomVOD


end module CustomVOD_Mod
