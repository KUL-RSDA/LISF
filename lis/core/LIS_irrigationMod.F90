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
module LIS_irrigationMod
!BOP
!
! !MODULE: LIS_irrigationMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  11 Nov 2012: Sujay Kumar; Initial implementation
!  29 May 2019; Jessica Erlingis; Incorporate Wanshu Nie's max/min GVF update
!  23 Feb 2022; Sara Modanesi; Incorporate Growing season to avoid a double
!  option (i.e., based on GVF and based on dyn LAI for Noah-MP.v.3.6)
!
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_mpiMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_historyMod, only : LIS_writevar_spread
  use LIS_fileIOMod, only : LIS_create_irrspread_filename, &
                            LIS_create_output_directory

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_irrigation_init
  public :: LIS_irrigation_run
  public :: LIS_irrigation_output
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_irrig_state      !data structure containing irrigation variables
!EOP  

  type, private :: irrig_type_dec
     real               :: outInterval
     character*100      :: models_used
     logical            :: stats_file_open
  end type irrig_type_dec

  type(irrig_type_dec),allocatable :: LIS_irrig_struc(:)

  type(ESMF_State),    allocatable :: LIS_irrig_state(:)

contains

!BOP
! 
! !ROUTINE: LIS_irrigation_init
! \label{LIS_irrigation_init}
! 
! !DESCRIPTION:
!
! Allocates memory for data structures used for reading 
! irrigation datasets. The irrigationdepth field is updated by the external
! files. The irrigation water equivalent fields are expected to be set
! by the model. 
! 
! !INTERFACE:
  subroutine LIS_irrigation_init

! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   use LIS_timeMgrMod,   only : LIS_registerAlarm, LIS_parseTimeString
!EOP
    integer       :: n
    integer       :: status
    integer       :: rc
    integer       :: ios, nid
    character*100 :: temp
    character*10  :: time
    character*1   :: nestid(2)
    logical       :: file_exists
! ___________________________________________________

 !- Read in Config file irrigation inputs:

  ! Read in type of irrigation scheme selected (spray,flood,drip):
    call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_type,&
         label="Irrigation scheme:",default="none",rc=rc)
    call LIS_verify(rc,&
         'Irrigation scheme: option not specified in the config file')

    if( LIS_rc%irrigation_type .ne. "none" ) then 

       write(LIS_logunit,*) "[INFO] Irrigation scheme selected:  ",&
                             trim(LIS_rc%irrigation_type)
 
       allocate(LIS_irrig_state(LIS_rc%nnest))
       allocate(LIS_irrig_struc(LIS_rc%nnest))

     ! Frequency with which irrigation field is written out:
       call ESMF_ConfigGetAttribute(LIS_config,time,&
            label="Irrigation output interval:",rc=rc)
       call LIS_verify(rc,"Irrigation output interval: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation output interval:  ",time

     ! Threshold for which irrigation is triggered, like for flood irrigation:
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_thresh,&
            label="Irrigation threshold:",rc=rc)
       call LIS_verify(rc,"Irrigation threshold: not defined")
       write(LIS_logunit,*) "[INFO] and irrigation threshold:  ",&
                             LIS_rc%irrigation_thresh

     ! SM Feb 2022 add double option for the start of growing season
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%growing_season,&
            label="Growing season:",rc=rc)
       call LIS_verify(rc,"Growing season: not defined")
       write(LIS_logunit,*) "and growing season:  ",&
                             LIS_rc%growing_season
     !SM Feb 2022 end changes   

     ! Parameters to control the GVF threshold based on the range of GVF
     ! (shdmax-shdmin) for which sprinkler irrigation is triggered:(WN)
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GVFparam1,&
            label="Irrigation GVF parameter 1:",rc=rc)
       call LIS_verify(rc,"Irrigation GVF parameter 1: not defined")
       write(LIS_logunit,*) "and irrigation GVF parameter 1:  ",&
                             LIS_rc%irrigation_GVFparam1

       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GVFparam2,&
            label="Irrigation GVF parameter 2:",rc=rc)
       call LIS_verify(rc,"Irrigation GVF parameter 2: not defined")
       write(LIS_logunit,*) "and irrigation GVF parameter 2:  ",&
                             LIS_rc%irrigation_GVFparam2

     ! Max. soil layer depth for irrigation to reach to (available for flood only):
       LIS_rc%irrigation_mxsoildpth = 1
       if( LIS_rc%irrigation_type == "Flood" ) then
          call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_mxsoildpth,&
               label="Irrigation max soil layer depth:", default=1, rc=rc)
          call LIS_verify(rc,"Irrigation max soil layer depth: not defined")
          write(LIS_logunit,*) "[INFO]and irrigation max soil depth:  ",&
                                LIS_rc%irrigation_mxsoildpth
       endif

     ! JE Remove irrigated water from groundwater
       LIS_rc%irrigation_GWabstraction = 0 ! Default is no
       ! Need to add model sanity check here to make sure model contains GW (?)
       call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_GWabstraction,&
            label="Groundwater abstraction for irrigation:",default=0,rc=rc)
       call LIS_verify(rc,"Groundwater abstraction for irrigation: not defined")
       write(LIS_logunit,*) "[INFO]and irrigation withdrawn from GW:  ",&
                             LIS_rc%irrigation_GWabstraction

!------Wanshu----irrigation scheduling based on DVEG On--------
     LIS_rc%irrigation_dveg  = 0 ! Default is no
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_dveg,&
            label="Irrigation scheduling based on dynamic vegetation:",default=0,rc=rc)
       call LIS_verify(rc,"Irrigation scheduling based on dynamic vegetation: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation scheduling based on dynamic vegetation:  ",&
                             LIS_rc%irrigation_dveg
!------------------------------------------------------

!------Wanshu---GW abstraction based on irrigation groundwater ratio data--------
     LIS_rc%irrigation_SourcePartition  = 0 ! Default is no
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_SourcePartition,&
            label="Irrigation source water partition:",default=0,rc=rc)
       call LIS_verify(rc,"Irrigation source water partition: not defined")
       write(LIS_logunit,*) "[INFO] Irrigation source water partition:  ",&
                             LIS_rc%irrigation_SourcePartition
!------------------------------------------------------

!-------------------Louise B---Output irrigation ensemble spread------------------
     LIS_rc%irrigation_outspread = 0 ! Default is no
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_outspread,&
            label="Irrigation output ensemble spread:",default=0,rc=rc)
     
     if ((LIS_rc%irrigation_outspread.ne.0).and.(LIS_rc%irrigation_type.ne."Sprinkler")) then
          write(LIS_logunit,*) "Irrigation output ensemble spread only possible with Sprinkler"
          call LIS_endrun()
     endif

     if (LIS_rc%irrigation_outspread.ne.0) then
          call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%irrigation_ensemstype,&
               label="Irrigation ensemble spread type:",default="std",rc=rc)
          call LIS_verify(rc,'Irrigation ensemble spread type: not defined (either "max-min" or "std")')
     endif

!------------------------------------------------------

     ! Register irrigation output interval:
       do n=1,LIS_rc%nnest
          call LIS_parseTimeString(time,LIS_irrig_struc(n)%outInterval)
          call LIS_registerAlarm("LIS irrigation output interval",&
               real(LIS_irrig_struc(n)%outInterval), &
               LIS_irrig_struc(n)%outInterval)
          LIS_irrig_struc(n)%models_used = trim(LIS_rc%irrigation_type)
          LIS_irrig_struc(n)%stats_file_open = .true.
       enddo

       do n=1,LIS_rc%nnest
          write(unit=temp,fmt='(i2.2)') n
          read(unit=temp,fmt='(2a1)') nestid
          
          LIS_irrig_state(n) = ESMF_StateCreate(name="LSM Irrigation State"//&
               nestid(1)//nestid(2), rc=status)
          call LIS_verify(status, &
               "ESMF_StateCreate failed in LIS_irrigation_init")
       enddo

    !- Initiate the irrigation scheme selected in lis.config file:
       call irrigationschemeinit(trim(LIS_rc%irrigation_type)//char(0),&
            LIS_irrig_state)

    endif

  end subroutine LIS_irrigation_init

!BOP
! 
! !ROUTINE: LIS_irrigation_run
! \label{LIS_irrigation_run}
! 
! !INTERFACE:
  subroutine LIS_irrigation_run(n)
! !USES: 
    implicit none

! !ARGUMENTS: 
    integer  :: n 

! !DESCRIPTION:
! This routine runs the specified irrigation model.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP

    if(LIS_rc%irrigation_type.ne."none") then     

       call getirrigationlsmstates(trim(LIS_rc%lsm)//"+"//&
            trim(LIS_rc%irrigation_type)//char(0), n,LIS_irrig_state(n))
       call applyirrigationupdates(trim(LIS_rc%irrigation_type)//char(0),&
            n,LIS_irrig_state(n))
       
    endif

  end subroutine LIS_irrigation_run

!BOP
! 
! !ROUTINE: LIS_irrigation_output
! \label{LIS_irrigation_output}
! 
! !INTERFACE:
  subroutine LIS_irrigation_output(n)
! !USES: 

    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_historyMod, only : LIS_writeModelOutput
    use LIS_fileIOMod,  only : LIS_create_output_directory, &
         LIS_create_output_filename,  &
         LIS_create_stats_filename
    use LIS_constantsMod, only : LIS_CONST_PATH_LEN

! !ARGUMENTS: 
    integer, intent(in)   :: n 

! !DESCRIPTION:
! This routine writes the irrigation model output.
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!     index of the nest
!  \end{description}
!EOP
    
    logical           :: alarmCheck,open_stats
    character(len=LIS_CONST_PATH_LEN) :: outfile, statsfile

    if(LIS_rc%irrigation_type.ne."none") then 
       alarmCheck = LIS_isAlarmRinging(LIS_rc,&
            "LIS irrigation output interval")
       if(alarmCheck) then 
          open_stats = .false. 
          if(LIS_rc%wopt.ne."none") then 
             if(LIS_masterproc) then 
                call LIS_create_output_directory('IRRIGATION')
                if (LIS_irrig_struc(n)%stats_file_open) then
                   call LIS_create_stats_filename(n,statsfile,"IRRIGATION")
                   LIS_irrig_struc(n)%stats_file_open = .false.
                   open_stats = .true.
                endif
             endif

             call LIS_create_output_filename(n,outfile,&
                  model_name ="IRRIGATION")

             call LIS_writeModelOutput(n,outfile,statsfile,              &
                  open_stats,outInterval=LIS_irrig_struc(n)%outInterval, &
                  nsoillayers=1, lyrthk = (/1.0/),                       &
                  nsoillayers2=1,                                        &
                  model_name=LIS_irrig_struc(n)%models_used,group=4)

             if(LIS_rc%irrigation_outspread.eq.1) then
               call writeEnsembleSpread_irrigation(n, LIS_irrig_state(n))
             endif

          endif
       endif
    endif
    
  end subroutine LIS_irrigation_output

!BOP
! 
! !ROUTINE: writeEnsembleSpread_irrigation
! \label{writeEnsembleSpread_irrigation}
!
! !INTERFACE: 
  subroutine writeEnsembleSpread_irrigation(n, irrigState)
!
! !DESCRIPTION: 
!  This routine writes the ensemble spread (standard deviation) 
!  of the irrigation amounts per day (kg m-2 d-1)
!
!  The arguments and variables are: 
!  \begin{description}
!   \item[n]    index of the nest 
!  \end{description}

!EOP
    integer,  intent(in)    :: n
    type(ESMF_State)        :: irrigState

    integer                :: ftn 
    integer                :: t
    character(len=LIS_CONST_PATH_LEN) :: spreadfile
    integer                :: shuffle, deflate, deflate_level
    integer                :: dimID(3)
    integer                :: irrigSumRate_id
    character*100          :: varname, vardimname, standard_name
    integer                :: status, rc, ierr
    type(ESMF_Field)       :: irrigSumRateField
    real,  pointer         :: irrigsumRate(:)


     ! Get irrigSumRate
     call ESMF_StateGet(irrigState, "Irrigation sum rate",irrigSumRateField,rc=rc)
     call LIS_verify(rc,'ESMF_StateGet failed for Irrigation sum rate in writeEnsembleSpread_irrigation')    
     call ESMF_FieldGet(irrigSumRateField, localDE=0,farrayPtr=irrigSumRate,rc=rc)
     call LIS_verify(rc,'ESMF_FieldGet failed for Irrigation sum rate in writeEnsembleSpread_irrigation')
          
     if(LIS_masterproc) then
          call LIS_create_output_directory('IRRIGATION')
          call LIS_create_irrspread_filename(n,spreadfile,&
               'IRRIGATION')

#if (defined USE_NETCDF4)
          status = nf90_create(path=spreadfile,cmode=nf90_hdf5,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in irrigation_Mod')
#endif
#if (defined USE_NETCDF3)
          status = nf90_create(path=spreadfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LIS_verify(status,&
               'creating netcdf file '//trim(spreadfile)//&
               ' failed in irrigation_Mod')
#endif

          if(LIS_rc%wopt.eq."1d gridspace") then 
               call LIS_verify(nf90_def_dim(ftn,'ngrid',&
                    LIS_rc%glbngrid_red(n),&
                    dimID(1)),'nf90_def_dim for ngrid failed in irrigation_mod')
          elseif(LIS_rc%wopt.eq."2d gridspace") then 
               call LIS_verify(nf90_def_dim(ftn,'east_west',LIS_rc%gnc(n),&
                    dimID(1)),'nf90_def_dim for east_west failed in irrigation_mod')
               call LIS_verify(nf90_def_dim(ftn,'north_south',LIS_rc%gnr(n),&
                    dimID(2)),'nf90_def_dim for north_south failed in irrigation_mod')
          endif

          call LIS_verify(nf90_put_att(ftn,&
               NF90_GLOBAL,"missing_value", LIS_rc%udef),&
               'nf90_put_att for missing_value failed in irrigation_mod')

!--------------------------------------------------------------------------
!  Ensemble spread -meta data
!--------------------------------------------------------------------------

          varname = "ensspread_IrrigationRate_daily"
          vardimname = "ensspread_IrrigationRate_daily"
          standard_name = "Ensemble_spread_for_IrrigationRate_daily"

          if(LIS_rc%wopt.eq."1d gridspace") then            
               call LIS_verify(nf90_def_var(ftn,varname,&
                    nf90_float,&
                    dimids = dimID(1), varID=IrrigSumRate_Id),&
                    'nf90_def_var for ensspread failed in irrigation_mod')
               
          elseif(LIS_rc%wopt.eq."2d gridspace") then 
               call LIS_verify(nf90_def_var(ftn,varname,&
                    nf90_float,&
                    dimids = dimID(1:2), varID=IrrigSumRate_Id),&
                    'nf90_def_var for ensspread failed in irrigation_mod')
          endif

#if(defined USE_NETCDF4)
          call LIS_verify(nf90_def_var_deflate(ftn,&
               IrrigSumRate_Id,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for ensspread failed in irrigation_mod')             
#endif
          call LIS_verify(nf90_put_att(ftn,IrrigSumRate_Id,&
               "standard_name",standard_name),&
               'nf90_put_att for ensspread failed in irrigation_mod')
          call LIS_verify(nf90_enddef(ftn),&
               'nf90_enddef failed in irrigation_mod')
     endif

     call LIS_writevar_spread(ftn,n,LIS_rc%lsm_index,IrrigSumRate_id, &
          irrigSumRate,1,LIS_rc%irrigation_ensemstype)

     if(LIS_masterproc) then 
          call LIS_verify(nf90_close(ftn),&
               'nf90_close failed in irrigation_mod')
     endif
     ! Reset sum
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
          irrigSumRate(t) = 0
     enddo


  end subroutine writeEnsembleSpread_irrigation

end module LIS_irrigationMod
