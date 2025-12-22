!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module AquaCrop_parmsMod
!BOP
!
! !MODULE: AquaCrop_parmsMod
!
! !DESCRIPTION:
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the AquaCrop
!  crop type from the AC_Crop.Inventory
!
! !REVISION HISTORY:
!
!  10 May 2024; Michel Becthold, Louise Busschaert, initial implementation
!
  use LDT_constantsMod, only: LDT_CONST_PATH_LEN
  use LDT_paramDataMod, only: LDT_paramEntry

  implicit none

  PRIVATE
  !------------------------------------------------------------------------------
  ! !PUBLIC MEMBER FUNCTIONS:
  !------------------------------------------------------------------------------
  public :: AquaCropParms_init    !allocates memory for required structures
  public :: AquaCropParms_writeHeader
  public :: AquaCropParms_writeData

  !------------------------------------------------------------------------------
  ! !PUBLIC TYPES:
  !------------------------------------------------------------------------------
  public :: AquaCrop_struc

  type, public :: aquacrop_type_dec
     ! -  AquaCrop LSM-specific:
     type(LDT_paramEntry) :: cropt   ! crop type
     type(LDT_paramEntry) :: comp_size ! compartment size
     type(LDT_paramEntry), allocatable :: tmin_cli(:) ! tmin climatology
     type(LDT_paramEntry), allocatable :: tmax_cli(:) ! tmax climatology
     integer :: nlayers !  number of soil layers
     real :: lthickness(5) ! thickness of layers, max 5 layers for AC
     integer :: max_comp ! fixed to 12
     integer :: tempcli_refyr ! reference year for cli, if nmetforcings>1, set to same value
     character(len=LDT_CONST_PATH_LEN), allocatable :: tempclimdir(:)
     character(len=LDT_CONST_PATH_LEN), allocatable :: tempclimfile(:)
     character(125), allocatable :: tempclim_gridtransform(:)
  end type aquacrop_type_dec

  type(aquacrop_type_dec), allocatable :: AquaCrop_struc(:)

contains

  subroutine AquaCropParms_init(flag)

    ! !USES:
    use ESMF
    use LDT_coreMod, only: LDT_rc, LDT_config
    use LDT_logMod,  only: LDT_verify, LDT_logunit, LDT_endrun
    !
    ! !DESCRIPTION:
    !
    ! Allocates memory for data structures for reading 
    ! the required AquaCrop datasets 
    ! 
    !  The routines invoked are: 
    !  \begin{description}
    !   \item[read_CONSTANT_AC_crop]\newline
    !    calls the reader of the crop type via the Crop.Inventory
    !    (for now: only CONSTANT method is implemented) 
    !   \item[define_AC_compartments]\newline
    !    calls the calculation fo the compartment number/size
    !   \item[read_AC_Tclim]\newline
    !    calls the reader of the Tmin and Tmax climatologies
    !  \end{description}
    !
    !EOP
    implicit none

    integer, intent(in) :: flag

    integer  :: n,k,m
    integer  :: rc

    external :: read_CONSTANT_AC_crop
    external :: define_AC_compartments
    external :: read_AC_Tclim

    character*3 :: months(12)
    data months /'jan','feb','mar','apr','may','jun','jul','aug', &
         'sep','oct','nov','dec'/

    ! _____________________________________________________________________

    allocate(AquaCrop_struc(LDT_rc%nnest))
    do n=1,LDT_rc%nnest
       AquaCrop_struc(n)%max_comp = 12
       ! Set crop type
       call set_param_attribs(Aquacrop_struc(n)%cropt, "AC_CROPT",&
            units="-", &
            full_name="Aquacrop crop type")
       call ESMF_ConfigFindLabel(LDT_config,"AquaCrop crop type data source:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,AquaCrop_struc(n)%cropt%source,rc=rc)
       call LDT_verify(rc,"AquaCrop crop type data source: not defined")
       allocate(Aquacrop_struc(n)%cropt%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            Aquacrop_struc(n)%cropt%vlevels))
       select case (AquaCrop_struc(n)%cropt%source)
       case( "CONSTANT" )
          call read_CONSTANT_AC_crop(&
               n,AquaCrop_struc(n)%cropt%value(:,:,1))
       case default
          write(LDT_logunit,*) "[ERR] crop type data source not valid for AquaCrop."
          write(LDT_logunit,*) "  Please select: CONSTANT"
          write(LDT_logunit,*) "Program stopping ..."
          call LDT_endrun
       end select
       ! End crop type

       ! Define compartment size
       call set_param_attribs(Aquacrop_struc(n)%comp_size, "AC_comp_size",&
            units="m", &
            full_name="Aquacrop compartment size")
       Aquacrop_struc(n)%comp_size%vlevels = AquaCrop_struc(n)%max_comp
       Aquacrop_struc(n)%comp_size%num_bins = AquaCrop_struc(n)%max_comp
       allocate(Aquacrop_struc(n)%comp_size%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            Aquacrop_struc(n)%comp_size%vlevels))
       call define_AC_compartments(n, AquaCrop_struc(n)%comp_size%value(:,:,:))


       !! Read temperature climatology file
       ! allocate variables with number of frocing sources
       allocate(AquaCrop_struc(n)%tempclimdir(LDT_rc%nmetforc))
       allocate(AquaCrop_struc(n)%tempclim_gridtransform(LDT_rc%nmetforc))
       allocate(Aquacrop_struc(n)%tmin_cli(LDT_rc%nmetforc))
       allocate(Aquacrop_struc(n)%tmax_cli(LDT_rc%nmetforc))

       ! Read options from ldt.config
       call ESMF_ConfigFindLabel(LDT_config,"AquaCrop temperature climatology directory:",rc=rc)
       call LDT_verify(rc,"AquaCrop temperature climatology directory: not defined")
       do m=1,LDT_rc%nmetforc
          call ESMF_ConfigGetAttribute(LDT_config,AquaCrop_struc(n)%tempclimdir(m),rc=rc)
       enddo

       call ESMF_ConfigFindLabel(LDT_config,"AquaCrop reference year for climatology:",rc=rc)
       call LDT_verify(rc,"AquaCrop reference year for climatology: not defined")
       call ESMF_ConfigGetAttribute(LDT_config,AquaCrop_struc(n)%tempcli_refyr,rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,"AquaCrop temperature climatology spatial transform:",rc=rc)
       call LDT_verify(rc,"AquaCrop temperature climatology spatial transform: not defined")
       do m=1,LDT_rc%nmetforc
          call ESMF_ConfigGetAttribute(LDT_config,AquaCrop_struc(n)%tempclim_gridtransform(m),rc=rc)
       enddo

       LDT_rc%monthlyData(n) = .true.

       do m = 1, LDT_rc%nmetforc
         ! tmin
         call set_param_attribs(Aquacrop_struc(n)%tmin_cli(m), "AC_Tmin_clim_"//trim(LDT_rc%metforc(m)),&
               units="K", &
               full_name="minimum temperature climatology ("//trim(LDT_rc%metforc(m))//")")
         Aquacrop_struc(n)%tmin_cli(m)%vlevels = 12
         Aquacrop_struc(n)%tmin_cli(m)%num_bins = 12
         allocate(Aquacrop_struc(n)%tmin_cli(m)%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               Aquacrop_struc(n)%tmin_cli(m)%vlevels))

         ! tmax
         call set_param_attribs(Aquacrop_struc(n)%tmax_cli(m), "AC_Tmax_clim_"//trim(LDT_rc%metforc(m)),&
               units="K", &
               full_name="maximum temperature climatology ("//trim(LDT_rc%metforc(m))//")")
         Aquacrop_struc(n)%tmax_cli(m)%vlevels = 12
         Aquacrop_struc(n)%tmax_cli(m)%num_bins = 12
         allocate(Aquacrop_struc(n)%tmax_cli(m)%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               Aquacrop_struc(n)%tmax_cli(m)%vlevels))

         ! Call function to read monthly data
         do k = 1,12
            AquaCrop_struc(n)%tempclimfile = &
                  trim(AquaCrop_struc(n)%tempclimdir(m))//'tmin.'//&
                  trim(months(k))//'.txt'
            call read_AC_Tclim(n, m, AquaCrop_struc(n)%tmin_cli(m)%value(:,:,k))
            AquaCrop_struc(n)%tempclimfile = &
                  trim(AquaCrop_struc(n)%tempclimdir(m))//'tmax.'//&
                  trim(months(k))//'.txt'
            call read_AC_Tclim(n, m, AquaCrop_struc(n)%tmax_cli(m)%value(:,:,k))
         enddo ! end months
      enddo ! end met forcing sources

    enddo ! End nest

  end subroutine AquaCropParms_init

  subroutine AquaCropParms_writeHeader(n,ftn,dimID,monthID)

    use LDT_historyMod, only: LDT_writeNetcdfDataHeader
    use LDT_logMod, only: LDT_verify
    use LDT_coreMod, only: LDT_rc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ftn
    integer, intent(in) :: dimID(3)
    integer, intent(in) :: monthID

    integer :: i, m
    integer   :: ndimID(3)  ! 3D, vlevel>1
    character(25) :: str

    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         Aquacrop_struc(n)%cropt)
    call LDT_verify(nf90_def_dim(ftn,'compartment',&
         AquaCrop_struc(n)%max_comp,ndimID(3)))
    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
         Aquacrop_struc(n)%comp_size)

    ndimID(3) = monthID
    do m = 1, LDT_rc%nmetforc
      call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
            Aquacrop_struc(n)%tmin_cli(m))
      call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
            Aquacrop_struc(n)%tmax_cli(m))
    enddo

    ! Add number of soil layers
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOIL_LAYERS", &
         AquaCrop_struc(n)%nlayers))
    ! Add thickness of soil layers
    do i=1,AquaCrop_struc(n)%nlayers
       write (str, '(i0)') i
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"THICKNESS_LAYER_"//trim(str),&
            AquaCrop_struc(n)%lthickness(i)))
    enddo

    ! Add Reference year for temperature climatology
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"AC_CLIM_REF_YEAR"), &
       AquaCrop_struc(n)%tempcli_refyr)

#endif

  end subroutine AquaCropParms_writeHeader

  subroutine AquaCropParms_writeData(n,ftn)

    use LDT_historyMod, only: LDT_writeNetcdfData
    use LDT_coreMod, only: LDT_rc

    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ftn

    integer :: m

    call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%cropt)
    call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%comp_size)

    do m = 1, LDT_rc%nmetforc
      call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%tmin_cli(m))
      call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%tmax_cli(m))
    enddo

  end subroutine AquaCropParms_writeData

  !BOP
  ! !ROUTINE:  set_param_attribs
  ! \label{set_param_attribs}
  !
  ! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, &
       units, full_name )

    ! !DESCRIPTION:
    !   This routine reads over the parameter attribute entries
    !   in the param_attribs.txt file.
    !
    ! !USES:
    use LDT_paramDataMod, only: LDT_paramEntry

    implicit none

    type(LDT_paramEntry),intent(inout) :: paramEntry
    character(len=*),    intent(in)    :: short_name
    character(len=*),     optional     :: units
    character(len=*),     optional     :: full_name

    character(20) :: unit_temp
    character(100):: name_temp
    ! ____________________________________________________

    if(present(units)) then
       unit_temp = units
    else
       unit_temp = "none"
    endif
    if(present(full_name)) then
       name_temp = full_name
    else
       name_temp = trim(short_name)
    endif

    paramEntry%short_name = trim(short_name)
    paramEntry%vlevels = 1
    paramEntry%selectOpt = 1
    paramEntry%source = "AquaCrop"
    paramEntry%units = unit_temp
    paramEntry%num_times = 1
    paramEntry%num_bins = 1
    paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module AquaCrop_parmsMod
