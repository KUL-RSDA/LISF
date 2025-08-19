! This defines the CustomLAI reader, based on the CustomNcReader
! It is basically just an instance of CustomNcReader_dec, and the associated
! subroutines call the respective ones of CustomNcReader_Mod.
! Normally it should be enough to copy, replace LAI by your new variable name,
! and potentially adapt the settings in CustomLAI_setup.
module CustomLAI_Mod
    use CustomNcReader_Mod, only: CustomNcReader_dec

    implicit none

    public :: CustomLAI_setup, read_CustomLAI, write_CustomLAI
    public :: CustomLAI_struc

    ! declare public reader array
    type(CustomNcReader_dec), allocatable :: CustomLAI_struc(:)

contains

    subroutine CustomLAI_setup(k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: CustomNcReader_setup

        implicit none

        ! !ARGUMENTS:
        integer                   :: k
        type(ESMF_State)          :: OBS_State(LIS_rc%nnest)
        type(ESMF_State)          :: OBS_Pert_State(LIS_rc%nnest)

        integer :: n

        allocate(CustomLAI_struc(LIS_rc%nnest))
        do n=1,LIS_rc%nnest
            CustomLAI_struc(n)%obsid = "Custom LAI"
            CustomLAI_struc(n)%varname = "LAI"
            CustomLAI_struc(n)%min_value = 0.0001
            CustomLAI_struc(n)%max_value = 10.0
            CustomLAI_struc(n)%qcmin_value = 0.0
            CustomLAI_struc(n)%qcmax_value = 100.0
        enddo

        call CustomNcReader_setup(CustomLAI_struc, k, OBS_State, OBS_Pert_State)

    end subroutine CustomLAI_setup

    subroutine read_CustomLAI(n, k, OBS_State, OBS_Pert_State)
        use ESMF, only: ESMF_State
        use LIS_coreMod, only: LIS_rc
        use CustomNcReader_Mod, only: read_CustomNetCDF
        integer, intent(in)       :: n, k
        type(ESMF_State)          :: OBS_State
        type(ESMF_State)          :: OBS_Pert_State
        call read_CustomNetCDF(CustomLAI_struc, n, k, OBS_State, OBS_Pert_State)
    end subroutine read_CustomLAI

    subroutine write_CustomLAI(n, k, OBS_State)
        use ESMF
        use CustomNcReader_Mod, only: write_CustomNetCDF
        integer,     intent(in)  :: n, k
        type(ESMF_State)         :: OBS_State
        call write_CustomNetCDF(n, k, OBS_State)
    end subroutine write_CustomLAI


end module CustomLAI_Mod
