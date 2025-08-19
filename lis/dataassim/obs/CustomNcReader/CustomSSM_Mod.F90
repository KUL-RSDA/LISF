! This defines the CustomSSM reader, based on the CustomNcReader
! It is basically just an instance of CustomNcReader_dec, and the associated
! subroutines call the respective ones of CustomNcReader_Mod.
! Normally it should be enough to copy, replace SSM by your new variable name,
! and potentially adapt the settings in CustomSSM_setup.
module CustomSSM_Mod
    use CustomNcReader_Mod, only: CustomNcReader_dec

    implicit none

    public :: CustomSSM_setup, read_CustomSSM, write_CustomSSM
    public :: CustomSSM_struc

    ! declare public reader array
    type(CustomNcReader_dec), allocatable :: CustomSSM_struc(:)

contains

    subroutine CustomSSM_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomSSM_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomSSM_struc(n)%obsid = "Custom SSM"
            CustomSSM_struc(n)%varname = "SSM"
            CustomSSM_struc(n)%min_value = 0.0001
            CustomSSM_struc(n)%max_value = 0.45
            CustomSSM_struc(n)%qcmin_value = 0.0
            CustomSSM_struc(n)%qcmax_value = 0.99
        enddo

        call CustomNcReader_setup(CustomSSM_struc, k, OBS_State, OBS_Pert_State)

    end subroutine CustomSSM_setup

    subroutine read_CustomSSM(n, k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: read_CustomNetCDF
        integer, intent(in)       :: n, k
        type(ESMF_State)          :: OBS_State
        type(ESMF_State)          :: OBS_Pert_State
        call read_CustomNetCDF(CustomSSM_struc, n, k, OBS_State, OBS_Pert_State)
    end subroutine read_CustomSSM

    subroutine write_CustomSSM(n, k, OBS_State)
        use ESMF
        use CustomNcReader_Mod, only: write_CustomNetCDF
        integer,     intent(in)  :: n, k
        type(ESMF_State)         :: OBS_State
        call write_CustomNetCDF(n, k, OBS_State)
    end subroutine write_CustomSSM


end module CustomSSM_Mod
