!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac70_main
! \label{Ac70_main}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   9/4/14: Shugong Wang; initial implementation for Ac70 with LIS-7
!   2/7/18: Soni Yatheendradas; code added for OPTUE to work
!
! !INTERFACE:
subroutine Ac70_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use Ac70_lsmMod
   !use other modules
    use ESMF
    use LIS_routingMod, only : LIS_runoff_state
    !!! MB_AC70
    use ac_global, only: typeproject_typeprm, &
                         typeproject_typepro, &
                         GetSimulParam_ThicknessTopSWC, &
                         GetRootZoneWC_Actual,&
                         GetRootZoneWC_FC,&
                         GetRootZoneWC_WP,&
                         GetRootZoneWC_SAT,&
                         GetRootZoneWC_Leaf,&
                         GetRootZoneWC_Thresh,&
                         GetRootZoneWC_Sen,&
                         GetRootZoneWC_ZtopAct,&
                         GetRootZoneWC_ZtopFC,&
                         GetRootZoneWC_ZtopWP,&
                         GetRootZoneWC_ZtopThresh,&
                         SetSimulParam_ThicknessTopSWC, &
                         SetRootZoneWC_Actual,&
                         SetRootZoneWC_FC,&
                         SetRootZoneWC_WP,&
                         SetRootZoneWC_SAT,&
                         SetRootZoneWC_Leaf,&
                         SetRootZoneWC_Thresh,&
                         SetRootZoneWC_Sen,&
                         SetRootZoneWC_ZtopAct,&
                         SetRootZoneWC_ZtopFC,&
                         SetRootZoneWC_ZtopWP,&
                         SetRootZoneWC_ZtopThresh,&
                         GetCompartment,&
                         SetCompartment,&
                         GetSoilLayer,&
                         SetSoilLayer,&
                         GetTotalSaltContent,&
                         GetTotalWaterContent,&
                         Geteffectiverain,&
                         GetSumWaBal,&
                         GetRootZoneSalt,&
                         GetSimulation,&
                         GetSimulation_SumGDD,&
                         GetSimulation_SumGDDfromDay1,&
                         SetTotalSaltContent,&
                         SetTotalWaterContent,&
                         Seteffectiverain,&
                         SetSumWaBal,&
                         SetRootZoneSalt,&
                         SetSimulation,&
                         GetIrrigation,&
                         SetIrrigation,&
                         GetCompartment_theta,&
                         SetCompartment_theta,&
                         GetIrriECw,&
                         GetManagement,&
                         GetPerennialPeriod,&
                         GetSimulParam,&
                         GetManagement_Cuttings,&
                         GetOnset,&
                         GetEndSeason,&
                         GetCrop,&
                         GetSoil,&
                         GetTemperatureRecord,&
                         GetClimRecord,&
                         GetRainRecord,&
                         GetEToRecord,&
                         SetIrriECw,&
                         SetManagement,&
                         SetPerennialPeriod,&
                         SetSimulParam,&
                         SetManagement_Cuttings,&
                         SetOnset,&
                         SetEndSeason,&
                         SetCrop,&
                         SetSoil,&
                         SetTemperatureRecord,&
                         SetClimRecord,&
                         SetRainRecord,&
                         SetEToRecord,&
                         GetSimulParam_GDDMethod,&

                         GetGenerateTimeMode,&
                        GetGenerateDepthMode,&
                        GetIrriMode,&
                        GetIrriMethod,&
                        GetDaySubmerged,&
                        GetMaxPlotNew,&
                        GetNrCompartments,&
                        GetIrriFirstDayNr,&
                        GetZiAqua,&
                        GetIniPercTAW,&
                        GetMaxPlotTr,&
                        GetOutputAggregate,&
                        GetEvapoEntireSoilSurface,&
                        GetPreDay,&
                        GetOutDaily,&
                        GetOut1Wabal,&
                        GetOut2Crop,&
                        GetOut3Prof,&
                        GetOut4Salt,&
                        GetOut5CompWC,&
                        GetOut6CompEC,&
                        GetOut7Clim,&
                        GetPart1Mult,&
                        GetPart2Eval,&
                        SetGenerateTimeMode,&
                        SetGenerateDepthMode,&
                        SetIrriMode,&
                        SetIrriMethod,&
                        SetDaySubmerged,&
                        SetMaxPlotNew,&
                        SetNrCompartments,&
                        SetIrriFirstDayNr,&
                        SetZiAqua,&
                        SetIniPercTAW,&
                        SetMaxPlotTr,&
                        SetOutputAggregate,&
                        SetEvapoEntireSoilSurface,&
                        SetPreDay,&
                        SetOutDaily,&
                        SetOut1Wabal,&
                        SetOut2Crop,&
                        SetOut3Prof,&
                        SetOut4Salt,&
                        SetOut5CompWC,&
                        SetOut6CompEC,&
                        SetOut7Clim,&
                        SetPart1Mult,&
                        SetPart2Eval,&


                        GetCCiActual,&
                        GetCCiprev,&
                        GetCCiTopEarlySen,&
                        GetCRsalt,&  
                        GetCRwater,& 
                        GetECdrain,& 
                        GetECiAqua,& 
                        GetECstorage,& 
                        GetEact,& 
                        GetEpot,& 
                        GetETo,&
                        GetDrain,&  
                        GetInfiltrated,&
                        GetRain,& 
                        GetRootingDepth,&
                        GetRunoff,& 
                        GetSaltInfiltr,&
                        GetSurf0,& 
                        GetSurfaceStorage,&
                        GetTact,&
                        GetTpot,&
                        GetTactWeedInfested,&
                        GetTmax,& 
                        GetTmin,& 
                        SetCCiActual,&
                        SetCCiprev,&
                        SetCCiTopEarlySen,&
                        SetCRsalt,&  
                        SetCRwater,& 
                        SetECdrain,& 
                        SetECiAqua,& 
                        SetECstorage,& 
                        SetEact,& 
                        SetEpot,& 
                        SetETo,&
                        SetDrain,&  
                        SetInfiltrated,&
                        SetRain,& 
                        SetRootingDepth,&
                        SetRunoff,& 
                        SetSaltInfiltr,&
                        SetSurf0,& 
                        SetSurfaceStorage,&
                        SetTact,&
                        SetTpot,&
                        SetTactWeedInfested,&
                        SetTmax,& 
                        SetTmin,&
                        GetIrriBeforeSeason,&
                        SetIrriBeforeSeason,&
                        GetIrriAfterSeason,&
                        SetIrriAfterSeason,&
                        GetCrop_Day1,&
                        DegreesDay,&
                         GetSimulation_ToDayNr, &
                         SetSimulation_ToDayNr, &
                        SetSimulation_SumGDDfromDay1,&
                        SetSimulation_SumGDD, &
                        SetClimRecord_DataType, &
                        SetClimRecord_fromd, &
                        SetClimRecord_fromdaynr, &
                        SetClimRecord_fromm, &
                        SetClimRecord_fromstring, &
                        SetClimRecord_fromy, &
                        SetClimRecord_NrObs, &
                        SetClimRecord_tod, &
                        SetClimRecord_todaynr, &
                        SetClimRecord_tom, &
                        SetClimRecord_tostring, &
                        SetClimRecord_toy,&
                        SetClimFile

           !!! MB_AC70
    !!! MB:
    use ac_project_input, only: ProjectInput 

    use ac_run, only:    SetDayNri,&
                         GetIrriInterval,&
                         GetIrriInfoRecord1,&
                         GetIrriInfoRecord2,&
                         SetIrriInterval,&
                         SetIrriInfoRecord1,&
                         SetIrriInfoRecord2,&
                         GetTheProjectFile,&
                        GetGwTable,&
                        GetPlotVarCrop,&
                        GetStressTot,&
                        GetCutInfoRecord1,&
                        GetCutInfoRecord2,&
                        GetTransfer,&
                        GetPreviousSum,&
                        GetTadj,&
                        GetGDDTadj,&
                        GetDayLastCut,&
                        GetNrCut,&
                        GetSumInterval,&
                        GetPreviousStressLevel,&
                        GetStressSFadjNEW,&
                        GetBin,&
                        GetBout,&
                        GetGDDayi,&
                        GetCO2i,&
                        GetFracBiomassPotSF,&
                        GetSumETo,&
                        GetSumGDD,&
                        GetZiprev,&
                        GetSumGDDPrev,&
                        GetCCxWitheredTpot,&
                        GetCCxWitheredTpotNoS,&
                        GetCoeffb0,&
                        GetCoeffb1,&
                        GetCoeffb2,&
                        GetCoeffb0Salt,&
                        GetCoeffb1Salt,&
                        GetCoeffb2Salt,&
                        GetStressLeaf,&
                        GetStressSenescence ,&
                        GetDayFraction,&
                        GetGDDayFraction,&
                        GetCGCref,&
                        GetGDDCGCref ,&
                        GetTimeSenescence ,&
                        GetSumKcTop,&
                        GetSumKcTopStress,&
                        GetSumKci,&
                        GetCCoTotal,&
                        GetCCxTotal,&
                        GetCDCTotal,&
                        GetGDDCDCTotal,&
                        GetCCxCropWeedsNoSFstress,&
                        GetWeedRCi,&
                        GetCCiActualWeedInfested,&
                        GetfWeedNoS,&
                        GetZeval,&
                        GetBprevSum,&
                        GetYprevSum,&
                        GetSumGDDcuts,&
                        GetHItimesBEF,&
                        GetScorAT1,&
                        GetScorAT2,&
                        GetHItimesAT1,&
                        GetHItimesAT2,&
                        GetHItimesAT,&
                        GetalfaHI,&
                        GetalfaHIAdj,&
                        GetNextSimFromDayNr ,&
                        GetDayNr1Eval,&
                        GetDayNrEval,&
                        GetLineNrEval,&
                        GetPreviousSumETo,&
                        GetPreviousSumGDD,&
                        GetPreviousBmob,&
                        GetPreviousBsto,&
                        GetStageCode,&
                        GetPreviousDayNr,&
                        GetNoYear,&
                        GetWaterTableInProfile,&
                        GetStartMode,&
                        GetNoMoreCrop,&
                        GetCGCadjustmentAfterCutting,&

                        SetGwTable,&
                        SetPlotVarCrop,&
                        SetStressTot,&
                        SetCutInfoRecord1,&
                        SetCutInfoRecord2,&
                        SetTransfer,&
                        SetPreviousSum,&
                        SetTadj,&
                        SetGDDTadj,&
                        SetDayLastCut,&
                        SetNrCut,&
                        SetSumInterval,&
                        SetPreviousStressLevel,&
                        SetStressSFadjNEW,&
                        SetBin,&
                        SetBout,&
                        SetGDDayi,&
                        SetCO2i,&
                        SetFracBiomassPotSF,&
                        SetSumETo,&
                        SetSumGDD,&
                        SetZiprev,&
                        SetSumGDDPrev,&
                        SetCCxWitheredTpot,&
                        SetCCxWitheredTpotNoS,&
                        SetCoeffb0,&
                        SetCoeffb1,&
                        SetCoeffb2,&
                        SetCoeffb0Salt,&
                        SetCoeffb1Salt,&
                        SetCoeffb2Salt,&
                        SetStressLeaf,&
                        SetStressSenescence ,&
                        SetDayFraction,&
                        SetGDDayFraction,&
                        SetCGCref,&
                        SetGDDCGCref ,&
                        SetTimeSenescence ,&
                        SetSumKcTop,&
                        SetSumKcTopStress,&
                        SetSumKci,&
                        SetCCoTotal,&
                        SetCCxTotal,&
                        SetCDCTotal,&
                        SetGDDCDCTotal,&
                        SetCCxCropWeedsNoSFstress,&
                        SetWeedRCi,&
                        SetCCiActualWeedInfested,&
                        SetfWeedNoS,&
                        SetZeval,&
                        SetBprevSum,&
                        SetYprevSum,&
                        SetSumGDDcuts,&
                        SetHItimesBEF,&
                        SetScorAT1,&
                        SetScorAT2,&
                        SetHItimesAT1,&
                        SetHItimesAT2,&
                        SetHItimesAT,&
                        SetalfaHI,&
                        SetalfaHIAdj,&
                        SetNextSimFromDayNr ,&
                        SetDayNr1Eval,&
                        SetDayNrEval,&
                        SetLineNrEval,&
                        SetPreviousSumETo,&
                        SetPreviousSumGDD,&
                        SetPreviousBmob,&
                        SetPreviousBsto,&
                        SetStageCode,&
                        SetPreviousDayNr,&
                        SetNoYear,&
                        SetWaterTableInProfile,&
                        SetStartMode,&
                        SetNoMoreCrop,&
                        SetCGCadjustmentAfterCutting,&
                        AdvanceOneTimeStep, &
                        ReadClimateNextDay, &
                        SetGDDVariablesNextDay, &
                         FinalizeRun1, &
                         FinalizeRun2, &
                         GetDayNri,&
                         GetCrop_Tbase, &
                         GetCrop_Tupper, &
                         FinalizeSimulation, &
                         InitializeSimulation, &
                         InitializeRunPart1, &
                         InitializeRunPart2, &
                         InitializeSimulationRunPart2, &
                         InitializeClimate
                         
    use ac_startunit, only:  FinalizeTheProgram, &
                         GetListProjectsFile, &
                         GetNumberOfProjects, &
                         GetProjectFileName, &
                         GetProjectType, &
                         GetSimulation_NrRuns, &
                         InitializeTheProgram, &
                         InitializeProject, &
                         WriteProjectsInfo


    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        dp,&
                        sp
    !!! MB_AC70

    implicit none

    !!! MB_AC70
    integer :: daynr, todaynr, iproject, nprojects
    logical :: ListProjectFileExist
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile

    !!! MB_AC70

! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    integer              :: itemp, countertemp
    real                 :: dt
    real                 :: lat, lon
    real                 :: Tmin_movmean
    real                 :: Tmin_mplr
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck

    integer               :: status
    integer               :: c,r,l
    integer               :: ios, nid,rivid,fldid

    integer            :: tid
!
! !DESCRIPTION:
!  This is the entry point for calling the Ac70 physics.
!  This routine calls the {\tt ac\_driver\_70 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for Ac70
    
    !MB: AC70
    real                 :: tmp_PREC_ac        ! 
    real                 :: tmp_TMIN_ac        ! 
    real                 :: tmp_TMAX_ac        ! 
    real                 :: tmp_ETo_ac        ! 


    ! check Ac70 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Ac70 model alarm")
    if (alarmCheck) Then

       do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          dt = LIS_rc%ts
          row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
          col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
          lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
          lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

          ! PREC_ac
          tmp_PREC_ac      = AC70_struc(n)%ac70(t)%PREC_ac  / AC70_struc(n)%forc_count
          ! TMIN_ac
          tmp_TMIN_ac      = AC70_struc(n)%ac70(t)%TMIN_ac  / AC70_struc(n)%forc_count
          ! TMAX_ac
          tmp_TMAX_ac      = AC70_struc(n)%ac70(t)%TMAX_ac  / AC70_struc(n)%forc_count
          ! ETo_ac
          tmp_ETo_ac      = AC70_struc(n)%ac70(t)%ETo_ac  / AC70_struc(n)%forc_count

          ! check validity of PREC_ac
          if(tmp_PREC_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable PREC_ac in Ac70"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            
          ! check validity of TMIN
          if(tmp_TMIN_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable TMIN in Ac70"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
          
          ! check validity of TMAX
          if(tmp_TMAX_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable TMAX in Ac70"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            
          ! check validity of ETo
          if(tmp_ETo_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable ETo in Ac70"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            

            !!! MB_AC70

            ! setting all global variables
            call SetRootZoneWC_Actual(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_Actual,8))
            call SetRootZoneWC_FC(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_FC,8))
            call SetRootZoneWC_WP(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_WP,8))
            call SetRootZoneWC_SAT(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_SAT,8))
            call SetRootZoneWC_Leaf(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_Leaf,8))
            call SetRootZoneWC_Thresh(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_Sen(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopAct(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopFC(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopWP(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopThresh(REAL(AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct,8))
            call SetCompartment(AC70_struc(n)%ac70(t)%Compartment)
            call SetTotalSaltContent(AC70_struc(n)%ac70(t)%TotalSaltContent)
            call SetTotalWaterContent(AC70_struc(n)%ac70(t)%TotalWaterContent)
            call Seteffectiverain(AC70_struc(n)%ac70(t)%effectiverain)
            call SetSumWaBal(AC70_struc(n)%ac70(t)%SumWaBal)
            call SetRootZoneSalt(AC70_struc(n)%ac70(t)%RootZoneSalt)
            call SetSimulation(AC70_struc(n)%ac70(t)%Simulation)
            call SetIrriInterval(AC70_struc(n)%ac70(t)%IrriInterval)
            call SetIrriInfoRecord1(AC70_struc(n)%ac70(t)%IrriInfoRecord1)
            call SetIrriInfoRecord2(AC70_struc(n)%ac70(t)%IrriInfoRecord2)
            call SetIrrigation(REAL(AC70_struc(n)%ac70(t)%Irrigation,8))
            do l=1, AC70_struc(n)%ac70(t)%NrCompartments
                 call SetCompartment_theta(l,REAL(AC70_struc(n)%ac70(t)%smc(l),8))
            enddo
            call SetIrriECw(AC70_struc(n)%ac70(t)%IrriECw) 
            call SetManagement(AC70_struc(n)%ac70(t)%Management) 
            call SetPerennialPeriod(AC70_struc(n)%ac70(t)%PerennialPeriod) 
            call SetSimulParam(AC70_struc(n)%ac70(t)%simulparam) 
            call SetManagement_Cuttings(AC70_struc(n)%ac70(t)%Cuttings) 
            call SetOnset(AC70_struc(n)%ac70(t)%onset) 
            call SetEndSeason(AC70_struc(n)%ac70(t)%endseason) 
            call SetCrop(AC70_struc(n)%ac70(t)%crop) 
            call SetSoil(AC70_struc(n)%ac70(t)%Soil) 
            call SetIrriBeforeSeason(AC70_struc(n)%ac70(t)%IrriBeforeSeason)
            call SetIrriAfterSeason(AC70_struc(n)%ac70(t)%IrriAfterSeason)
            call SetSoilLayer(AC70_struc(n)%ac70(t)%soillayer)
            call SetDayNri(AC70_struc(n)%ac70(t)%daynri)

            call SetGenerateTimeMode(AC70_struc(n)%ac70(t)%GenerateTimeMode) 
            call SetGenerateDepthMode(AC70_struc(n)%ac70(t)%GenerateDepthMode) 
            call SetIrriMode(AC70_struc(n)%ac70(t)%IrriMode) 
            call SetDaySubmerged(AC70_struc(n)%ac70(t)%DaySubmerged) 
            call SetMaxPlotNew(AC70_struc(n)%ac70(t)%MaxPlotNew) 
            call SetNrCompartments(AC70_struc(n)%ac70(t)%NrCompartments) 
            call SetIrriFirstDayNr(AC70_struc(n)%ac70(t)%IrriFirstDayNr) 
            call SetZiAqua(AC70_struc(n)%ac70(t)%ZiAqua) 
            call SetIniPercTAW(AC70_struc(n)%ac70(t)%IniPercTAW) 
            call SetMaxPlotTr(AC70_struc(n)%ac70(t)%MaxPlotTr)
            call SetOutputAggregate(AC70_struc(n)%ac70(t)%OutputAggregate) 

            call SetEvapoEntireSoilSurface(AC70_struc(n)%ac70(t)%EvapoEntireSoilSurface) 
            call SetPreDay(AC70_struc(n)%ac70(t)%PreDay) 
            call SetOutDaily(AC70_struc(n)%ac70(t)%OutDaily) 
            call SetOut1Wabal(AC70_struc(n)%ac70(t)%Out1Wabal) 
            call SetOut2Crop(AC70_struc(n)%ac70(t)%Out2Crop) 
            call SetOut3Prof(AC70_struc(n)%ac70(t)%Out3Prof) 
            call SetOut4Salt(AC70_struc(n)%ac70(t)%Out4Salt) 
            call SetOut5CompWC(AC70_struc(n)%ac70(t)%Out5CompWC) 
            call SetOut6CompEC(AC70_struc(n)%ac70(t)%Out6CompEC) 
            call SetOut7Clim(AC70_struc(n)%ac70(t)%Out7Clim) 
            call SetPart1Mult(AC70_struc(n)%ac70(t)%Part1Mult) 
            call SetPart2Eval(AC70_struc(n)%ac70(t)%Part2Eval) 

            !
            call SetCCiActual(AC70_struc(n)%ac70(t)%CCiActual)
            call SetCCiprev(AC70_struc(n)%ac70(t)%CCiprev)

            call SetCCiTopEarlySen(AC70_struc(n)%ac70(t)%CCiTopEarlySen)
            call SetCRsalt(AC70_struc(n)%ac70(t)%CRsalt)
            call SetCRwater(AC70_struc(n)%ac70(t)%CRwater)
            call SetECdrain(AC70_struc(n)%ac70(t)%ECdrain)
            call SetEciAqua(AC70_struc(n)%ac70(t)%ECiAqua)
            call SetECstorage(AC70_struc(n)%ac70(t)%ECstorage)
            call SetEact(AC70_struc(n)%ac70(t)%Eact)
            call SetEpot(AC70_struc(n)%ac70(t)%Epot)
            call SetDrain(AC70_struc(n)%ac70(t)%Drain)
            call SetInfiltrated(AC70_struc(n)%ac70(t)%Infiltrated)
            call SetRootingDepth(AC70_struc(n)%ac70(t)%RootingDepth)
            call SetRunoff(AC70_struc(n)%ac70(t)%Runoff)
            call SetSaltInfiltr(AC70_struc(n)%ac70(t)%SaltInfiltr)
            call SetSurf0(AC70_struc(n)%ac70(t)%Surf0)
            call SetSurfaceStorage(AC70_struc(n)%ac70(t)%SurfaceStorage)
            call SetTact(AC70_struc(n)%ac70(t)%Tact)
            call SetTpot(AC70_struc(n)%ac70(t)%Tpot)
            call SetTactWeedInfested(AC70_struc(n)%ac70(t)%TactWeedInfested)

            call SetGwTable(AC70_struc(n)%ac70(t)%GwTable)
            call SetPlotVarCrop(AC70_struc(n)%ac70(t)%PlotVarCrop)
            call SetStressTot(AC70_struc(n)%ac70(t)%StressTot)
            call SetCutInfoRecord1(AC70_struc(n)%ac70(t)%CutInfoRecord1)
            call SetCutInfoRecord2(AC70_struc(n)%ac70(t)%CutInfoRecord2)
            call SetTransfer(AC70_struc(n)%ac70(t)%Transfer)
            call SetPreviousSum(AC70_struc(n)%ac70(t)%PreviousSum)
            call SetTadj(AC70_struc(n)%ac70(t)%Tadj)
            call SetGDDTadj(AC70_struc(n)%ac70(t)%GDDTadj)
            call SetDayLastCut(AC70_struc(n)%ac70(t)%DayLastCut)
            call SetNrCut(AC70_struc(n)%ac70(t)%NrCut)
            call SetSumInterval(AC70_struc(n)%ac70(t)%SumInterval)
            call SetPreviousStressLevel(int(AC70_struc(n)%ac70(t)%PreviousStressLevel,kind=int32))
            call SetStressSFadjNEW(int(AC70_struc(n)%ac70(t)%StressSFadjNEW,kind=int32))
            call SetBin(AC70_struc(n)%ac70(t)%Bin)
            call SetBout(AC70_struc(n)%ac70(t)%Bout)
            call SetCO2i(AC70_struc(n)%ac70(t)%CO2i)
            call SetFracBiomassPotSF(AC70_struc(n)%ac70(t)%FracBiomassPotSF)
            call SetSumETo(AC70_struc(n)%ac70(t)%SumETo)
            call SetSumGDD(AC70_struc(n)%ac70(t)%SumGDD)
            call SetZiprev(AC70_struc(n)%ac70(t)%Ziprev)
            call SetSumGDDPrev(AC70_struc(n)%ac70(t)%SumGDDPrev)
            call SetCCxWitheredTpot(AC70_struc(n)%ac70(t)%CCxWitheredTpot)
            call SetCCxWitheredTpotNoS(AC70_struc(n)%ac70(t)%CCxWitheredTpotNoS)
            call SetCoeffb0(AC70_struc(n)%ac70(t)%Coeffb0)
            call SetCoeffb1(AC70_struc(n)%ac70(t)%Coeffb1)
            call SetCoeffb2(AC70_struc(n)%ac70(t)%Coeffb2)
            call SetCoeffb0Salt(AC70_struc(n)%ac70(t)%Coeffb0Salt)
            call SetCoeffb1Salt(AC70_struc(n)%ac70(t)%Coeffb1Salt)
            call SetCoeffb2Salt(AC70_struc(n)%ac70(t)%Coeffb2Salt)
            call SetStressLeaf(AC70_struc(n)%ac70(t)%StressLeaf)
            call SetStressSenescence(AC70_struc(n)%ac70(t)%StressSenescence)
            call SetDayFraction(AC70_struc(n)%ac70(t)%DayFraction)
            call SetGDDayFraction(AC70_struc(n)%ac70(t)%GDDayFraction)
            call SetCGCref(AC70_struc(n)%ac70(t)%CGCref)
            call SetGDDCGCref(AC70_struc(n)%ac70(t)%GDDCGCref)
            call SetTimeSenescence(AC70_struc(n)%ac70(t)%TimeSenescence)
            call SetSumKcTop(AC70_struc(n)%ac70(t)%SumKcTop)
            call SetSumKcTopStress(AC70_struc(n)%ac70(t)%SumKcTopStress)
            call SetSumKci(AC70_struc(n)%ac70(t)%SumKci)
            call SetCCoTotal(AC70_struc(n)%ac70(t)%CCoTotal)
            call SetCCxTotal(AC70_struc(n)%ac70(t)%CCxTotal)
            call SetCDCTotal(AC70_struc(n)%ac70(t)%CDCTotal)
            call SetGDDCDCTotal(AC70_struc(n)%ac70(t)%GDDCDCTotal)
            call SetCCxCropWeedsNoSFstress(AC70_struc(n)%ac70(t)%CCxCropWeedsNoSFstress)
            call SetWeedRCi(AC70_struc(n)%ac70(t)%WeedRCi)
            call SetCCiActualWeedInfested(AC70_struc(n)%ac70(t)%CCiActualWeedInfested)
            call SetfWeedNoS(AC70_struc(n)%ac70(t)%fWeedNoS)
            call SetZeval(AC70_struc(n)%ac70(t)%Zeval)
            call SetBprevSum(AC70_struc(n)%ac70(t)%BprevSum)
            call SetYprevSum(AC70_struc(n)%ac70(t)%YprevSum)
            call SetSumGDDcuts(AC70_struc(n)%ac70(t)%SumGDDcuts)
            call SetHItimesBEF(AC70_struc(n)%ac70(t)%HItimesBEF)
            call SetScorAT1(AC70_struc(n)%ac70(t)%ScorAT1)
            call SetScorAT2(AC70_struc(n)%ac70(t)%ScorAT2)
            call SetHItimesAT1(AC70_struc(n)%ac70(t)%HItimesAT1)
            call SetHItimesAT2(AC70_struc(n)%ac70(t)%HItimesAT2)
            call SetHItimesAT(AC70_struc(n)%ac70(t)%HItimesAT)
            call SetalfaHI(AC70_struc(n)%ac70(t)%alfaHI)
            call SetalfaHIAdj(AC70_struc(n)%ac70(t)%alfaHIAdj)
            call SetNextSimFromDayNr(AC70_struc(n)%ac70(t)%NextSimFromDayNr)
            call SetDayNr1Eval(AC70_struc(n)%ac70(t)%DayNr1Eval)
            call SetDayNrEval(AC70_struc(n)%ac70(t)%DayNrEval)
            call SetLineNrEval(int(AC70_struc(n)%ac70(t)%LineNrEval,kind=int32))
            call SetPreviousSumETo(AC70_struc(n)%ac70(t)%PreviousSumETo)
            call SetPreviousSumGDD(AC70_struc(n)%ac70(t)%PreviousSumGDD)
            call SetPreviousBmob(AC70_struc(n)%ac70(t)%PreviousBmob)
            call SetPreviousBsto(AC70_struc(n)%ac70(t)%PreviousBsto)
            call SetStageCode(AC70_struc(n)%ac70(t)%StageCode)
            call SetPreviousDayNr(AC70_struc(n)%ac70(t)%PreviousDayNr)
            call SetNoYear(AC70_struc(n)%ac70(t)%NoYear)
            call SetWaterTableInProfile(AC70_struc(n)%ac70(t)%WaterTableInProfile)
            call SetStartMode(AC70_struc(n)%ac70(t)%StartMode)
            call SetNoMoreCrop(AC70_struc(n)%ac70(t)%NoMoreCrop)
            call SetCGCadjustmentAfterCutting(AC70_struc(n)%ac70(t)%CGCadjustmentAfterCutting)
            !call SetSimulation_ToDayNr(AC70_struc(n)%ac70(t)%Simulation%ToDayNr)
            call SetGDDayi(AC70_struc(n)%ac70(t)%GDDayi)

            !!! initialize run (year)

             if (AC70_struc(n)%ac70(t)%InitializeRun .eq. 1) then
                ! Replaces LoadSimulationRunProject in LoadSimulationRunProject
                !call SetSimulation_YearSeason(YearSeason_temp)
                !call SetCrop_Day1(TempInt)
                !call SetCrop_DayN(TempInt)
                !SetCO2FileFull

                ! Run InitializeRunPart1 ( in future without LoadSimulationRunProject)
                        call SetClimRecord_DataType(0_int8)
                        call SetClimRecord_fromd(0)
                        call SetClimRecord_fromdaynr(ProjectInput(1)%Simulation_DayNr1)
                        call SetClimRecord_fromm(0)
                        call SetClimRecord_fromstring("")
                        call SetClimRecord_fromy(LIS_rc%syr)
                        call SetClimRecord_NrObs(999)
                        call SetClimRecord_tod(0)
                        call SetClimRecord_todaynr(ProjectInput(GetSimulation_NrRuns())%Simulation_DayNrN)
                        call SetClimRecord_tom(0)
                        call SetClimRecord_tostring("")
                        call SetClimRecord_toy(0)
                        call SetClimFile('(External)')
                        !call SetClimFile('EToRainTempFile')
                        !call SetClimDescription('Read ETo/RAIN/TEMP data set')

                AC70_struc(n)%ac70(t)%WPi = 0._dp

                call InitializeRunPart1(AC70_struc(n)%ac70(t)%irun, AC70_struc(n)%ac70(t)%TheProjectType);
                if (trim(LIS_rc%metforc(1)) == 'MERRA2_AC') then
                  call SetRain(real(tmp_PREC_ac,kind=dp))
                  call SetTmin(real(tmp_TMIN_ac,kind=dp))
                  call SetTmax(real(tmp_TMAX_ac,kind=dp))
                  call SetETo(real(tmp_ETo_ac,kind=dp))
                else ! read from AC input
                  call InitializeClimate();
                end if

                !call InitializeRunPart2(AC70_struc(n)%ac70(t)%irun, AC70_struc(n)%ac70(t)%TheProjectType);
                call InitializeSimulationRunPart2()
                AC70_struc(n)%ac70(t)%InitializeRun = 0
            end if

            ! for MERRA2_AC --> first set climate variables then advanceonetimestep
            call SetRain(real(tmp_PREC_ac,kind=dp))
            call SetTmin(real(tmp_TMIN_ac,kind=dp))
            call SetTmax(real(tmp_TMAX_ac,kind=dp))
            call SetETo(real(tmp_ETo_ac,kind=dp))
            ! Sum of GDD at end of first day
            call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                 GetTmax(), GetSimulParam_GDDMethod()))
            if (GetDayNri() >= GetCrop_Day1()) then
               if (GetDayNri() == GetCrop_Day1()) then
                  call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
               end if
               call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                  GetGDDayi())
            end if
            !write(LIS_logunit,*) &
            !                '[INFO] AdvanceOneTimeStep AquaCrop, day in month: ', LIS_rc
            call AdvanceOneTimeStep(AC70_struc(n)%ac70(t)%WPi)

            !AC70_struc(n)%ac70(t)%smc(1)= GetRootZoneWC_ZtopAct()/(GetSimulParam_ThicknessTopSWC()*10._dp)
            AC70_struc(n)%ac70(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
            AC70_struc(n)%ac70(t)%RootZoneWC_FC = GetRootZoneWC_FC()
            AC70_struc(n)%ac70(t)%RootZoneWC_WP = GetRootZoneWC_WP()
            AC70_struc(n)%ac70(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
            AC70_struc(n)%ac70(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
            AC70_struc(n)%ac70(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
            AC70_struc(n)%ac70(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
            AC70_struc(n)%ac70(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
            AC70_struc(n)%ac70(t)%Compartment = GetCompartment()
            AC70_struc(n)%ac70(t)%TotalSaltContent = GetTotalSaltContent()
            AC70_struc(n)%ac70(t)%TotalWaterContent = GetTotalWaterContent()
            AC70_struc(n)%ac70(t)%effectiverain = Geteffectiverain()
            AC70_struc(n)%ac70(t)%SumWaBal = GetSumWaBal()
            AC70_struc(n)%ac70(t)%RootZoneSalt = GetRootZoneSalt()
            AC70_struc(n)%ac70(t)%Simulation = GetSimulation()
            AC70_struc(n)%ac70(t)%IrriInterval = GetIrriInterval()
            AC70_struc(n)%ac70(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC70_struc(n)%ac70(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC70_struc(n)%ac70(t)%Irrigation = GetIrrigation()
            AC70_struc(n)%ac70(t)%IrriBeforeSeason = GetIrriBeforeSeason()
            AC70_struc(n)%ac70(t)%IrriAfterSeason = GetIrriAfterSeason()
            AC70_struc(n)%ac70(t)%SoilLayer = GetSoilLayer()
            AC70_struc(n)%ac70(t)%daynri = GetDayNri()
            do l=1, AC70_struc(n)%ac70(t)%NrCompartments
                 AC70_struc(n)%ac70(t)%smc(l) = GetCompartment_theta(l)
            enddo
            !write(*,'(e23.15e3)') AC70_struc(n)%ac70(t)%ac70smc(1)
            AC70_struc(n)%ac70(t)%IrriECw = GetIrriECw()
            AC70_struc(n)%ac70(t)%Management = GetManagement()
            AC70_struc(n)%ac70(t)%PerennialPeriod = GetPerennialPeriod()
            AC70_struc(n)%ac70(t)%simulparam = GetSimulParam()
            AC70_struc(n)%ac70(t)%Cuttings = GetManagement_Cuttings()
            AC70_struc(n)%ac70(t)%onset = GetOnset()
            AC70_struc(n)%ac70(t)%endseason = GetEndSeason()
            AC70_struc(n)%ac70(t)%crop = GetCrop()
            AC70_struc(n)%ac70(t)%Soil = GetSoil()
            AC70_struc(n)%ac70(t)%TemperatureRecord = GetTemperatureRecord()
            AC70_struc(n)%ac70(t)%ClimRecord = GetClimRecord()
            AC70_struc(n)%ac70(t)%RainRecord = GetRainRecord()
            AC70_struc(n)%ac70(t)%EToRecord = GetEToRecord()

            AC70_struc(n)%ac70(t)%GenerateTimeMode = GetGenerateTimeMode()
            AC70_struc(n)%ac70(t)%GenerateDepthMode = GetGenerateDepthMode()
            AC70_struc(n)%ac70(t)%IrriMode = GetIrriMode()
            AC70_struc(n)%ac70(t)%IrriMethod = GetIrriMethod()
            AC70_struc(n)%ac70(t)%DaySubmerged = GetDaySubmerged()
            AC70_struc(n)%ac70(t)%MaxPlotNew = GetMaxPlotNew()
            AC70_struc(n)%ac70(t)%NrCompartments = GetNrCompartments()
            AC70_struc(n)%ac70(t)%IrriFirstDayNr = GetIrriFirstDayNr()
            AC70_struc(n)%ac70(t)%ZiAqua = GetZiAqua()
            AC70_struc(n)%ac70(t)%IniPercTAW = GetIniPercTAW()
            AC70_struc(n)%ac70(t)%MaxPlotTr = GetMaxPlotTr()
            AC70_struc(n)%ac70(t)%OutputAggregate = GetOutputAggregate()

            AC70_struc(n)%ac70(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
            AC70_struc(n)%ac70(t)%PreDay = GetPreDay()
            AC70_struc(n)%ac70(t)%OutDaily = GetOutDaily()
            AC70_struc(n)%ac70(t)%Out1Wabal = GetOut1Wabal()
            AC70_struc(n)%ac70(t)%Out2Crop = GetOut2Crop()
            AC70_struc(n)%ac70(t)%Out3Prof = GetOut3Prof()
            AC70_struc(n)%ac70(t)%Out4Salt = GetOut4Salt()
            AC70_struc(n)%ac70(t)%Out5CompWC = GetOut5CompWC()
            AC70_struc(n)%ac70(t)%Out6CompEC = GetOut6CompEC()
            AC70_struc(n)%ac70(t)%Out7Clim = GetOut7Clim()
            AC70_struc(n)%ac70(t)%Part1Mult = GetPart1Mult()
            AC70_struc(n)%ac70(t)%Part2Eval = GetPart2Eval()

            !
            AC70_struc(n)%ac70(t)%CCiActual = GetCCiActual()
            AC70_struc(n)%ac70(t)%CCiprev = GetCCiprev()
            AC70_struc(n)%ac70(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC70_struc(n)%ac70(t)%CRsalt = GetCRsalt () ! gram/m2
            AC70_struc(n)%ac70(t)%CRwater = GetCRwater() ! mm/day
            AC70_struc(n)%ac70(t)%ECdrain = GetECdrain() ! EC drain water dS/m
            AC70_struc(n)%ac70(t)%ECiAqua = GetECiAqua() ! EC of the groundwater table in dS/m
            AC70_struc(n)%ac70(t)%ECstorage = GetECstorage() !EC surface storage dS/m
            AC70_struc(n)%ac70(t)%Eact = GetEact() ! mm/day
            AC70_struc(n)%ac70(t)%Epot = GetEpot() ! mm/day
            AC70_struc(n)%ac70(t)%ETo_ac = GetETo() ! mm/day
            AC70_struc(n)%ac70(t)%Drain = GetDrain()  ! mm/day
            AC70_struc(n)%ac70(t)%Infiltrated = GetInfiltrated() ! mm/day
            AC70_struc(n)%ac70(t)%PREC_ac = GetRain()  ! mm/day
            AC70_struc(n)%ac70(t)%RootingDepth = GetRootingDepth()
            AC70_struc(n)%ac70(t)%Runoff = GetRunoff()  ! mm/day
            AC70_struc(n)%ac70(t)%SaltInfiltr = GetSaltInfiltr() ! salt infiltrated in soil profile Mg/ha
            AC70_struc(n)%ac70(t)%Surf0 = GetSurf0()  ! surface water [mm] begin day
            AC70_struc(n)%ac70(t)%SurfaceStorage = GetSurfaceStorage() !mm/day
            AC70_struc(n)%ac70(t)%Tact = GetTact() ! mm/day
            AC70_struc(n)%ac70(t)%Tpot = GetTpot() ! mm/day
            AC70_struc(n)%ac70(t)%TactWeedInfested = GetTactWeedInfested() !mm/day
            AC70_struc(n)%ac70(t)%Tmax_ac = GetTmax() ! degC
            AC70_struc(n)%ac70(t)%Tmin_ac =GetTmin() ! degC


            AC70_struc(n)%ac70(t)%GwTable = GetGwTable()
            AC70_struc(n)%ac70(t)%PlotVarCrop = GetPlotVarCrop()
            AC70_struc(n)%ac70(t)%StressTot = GetStressTot()
            AC70_struc(n)%ac70(t)%CutInfoRecord1 = GetCutInfoRecord1()
            AC70_struc(n)%ac70(t)%CutInfoRecord2 = GetCutInfoRecord2()
            AC70_struc(n)%ac70(t)%Transfer = GetTransfer()
            AC70_struc(n)%ac70(t)%PreviousSum = GetPreviousSum()
            AC70_struc(n)%ac70(t)%Tadj = GetTadj()
            AC70_struc(n)%ac70(t)%GDDTadj = GetGDDTadj()
            AC70_struc(n)%ac70(t)%DayLastCut = GetDayLastCut()
            AC70_struc(n)%ac70(t)%NrCut = GetNrCut()
            AC70_struc(n)%ac70(t)%SumInterval = GetSumInterval()
            AC70_struc(n)%ac70(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC70_struc(n)%ac70(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC70_struc(n)%ac70(t)%Bin = GetBin()
            AC70_struc(n)%ac70(t)%Bout = GetBout()
            AC70_struc(n)%ac70(t)%GDDayi = GetGDDayi()
            AC70_struc(n)%ac70(t)%CO2i = GetCO2i()
            AC70_struc(n)%ac70(t)%FracBiomassPotSF = GetFracBiomassPotSF()
            AC70_struc(n)%ac70(t)%SumETo = GetSumETo()
            AC70_struc(n)%ac70(t)%SumGDD = GetSumGDD()
            AC70_struc(n)%ac70(t)%Ziprev = GetZiprev()
            AC70_struc(n)%ac70(t)%SumGDDPrev = GetSumGDDPrev()
            AC70_struc(n)%ac70(t)%CCxWitheredTpot = GetCCxWitheredTpot()
            AC70_struc(n)%ac70(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC70_struc(n)%ac70(t)%Coeffb0 = GetCoeffb0()
            AC70_struc(n)%ac70(t)%Coeffb1 = GetCoeffb1()
            AC70_struc(n)%ac70(t)%Coeffb2 = GetCoeffb2()
            AC70_struc(n)%ac70(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC70_struc(n)%ac70(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC70_struc(n)%ac70(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC70_struc(n)%ac70(t)%StressLeaf = GetStressLeaf()
            AC70_struc(n)%ac70(t)%StressSenescence = GetStressSenescence ()
            AC70_struc(n)%ac70(t)%DayFraction = GetDayFraction()
            AC70_struc(n)%ac70(t)%GDDayFraction = GetGDDayFraction()
            AC70_struc(n)%ac70(t)%CGCref = GetCGCref()
            AC70_struc(n)%ac70(t)%GDDCGCref = GetGDDCGCref ()
            AC70_struc(n)%ac70(t)%TimeSenescence = GetTimeSenescence ()
            AC70_struc(n)%ac70(t)%SumKcTop = GetSumKcTop()
            AC70_struc(n)%ac70(t)%SumKcTopStress = GetSumKcTopStress()
            AC70_struc(n)%ac70(t)%SumKci = GetSumKci()
            AC70_struc(n)%ac70(t)%CCoTotal = GetCCoTotal()
            AC70_struc(n)%ac70(t)%CCxTotal = GetCCxTotal()
            AC70_struc(n)%ac70(t)%CDCTotal = GetCDCTotal()
            AC70_struc(n)%ac70(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC70_struc(n)%ac70(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC70_struc(n)%ac70(t)%WeedRCi = GetWeedRCi()
            AC70_struc(n)%ac70(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC70_struc(n)%ac70(t)%fWeedNoS = GetfWeedNoS()
            AC70_struc(n)%ac70(t)%Zeval = GetZeval()
            AC70_struc(n)%ac70(t)%BprevSum = GetBprevSum()
            AC70_struc(n)%ac70(t)%YprevSum = GetYprevSum()
            AC70_struc(n)%ac70(t)%SumGDDcuts = GetSumGDDcuts()
            AC70_struc(n)%ac70(t)%HItimesBEF = GetHItimesBEF()
            AC70_struc(n)%ac70(t)%ScorAT1 = GetScorAT1()
            AC70_struc(n)%ac70(t)%ScorAT2 = GetScorAT2()
            AC70_struc(n)%ac70(t)%HItimesAT1 = GetHItimesAT1()
            AC70_struc(n)%ac70(t)%HItimesAT2 = GetHItimesAT2()
            AC70_struc(n)%ac70(t)%HItimesAT = GetHItimesAT()
            AC70_struc(n)%ac70(t)%alfaHI = GetalfaHI()
            AC70_struc(n)%ac70(t)%alfaHIAdj = GetalfaHIAdj()
            AC70_struc(n)%ac70(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
            AC70_struc(n)%ac70(t)%DayNr1Eval = GetDayNr1Eval()
            AC70_struc(n)%ac70(t)%DayNrEval = GetDayNrEval()
            AC70_struc(n)%ac70(t)%LineNrEval = GetLineNrEval()
            AC70_struc(n)%ac70(t)%PreviousSumETo = GetPreviousSumETo()
            AC70_struc(n)%ac70(t)%PreviousSumGDD = GetPreviousSumGDD()
            AC70_struc(n)%ac70(t)%PreviousBmob = GetPreviousBmob()
            AC70_struc(n)%ac70(t)%PreviousBsto = GetPreviousBsto()
            AC70_struc(n)%ac70(t)%StageCode = GetStageCode()
            AC70_struc(n)%ac70(t)%PreviousDayNr = GetPreviousDayNr()
            AC70_struc(n)%ac70(t)%NoYear = GetNoYear()
            AC70_struc(n)%ac70(t)%WaterTableInProfile = GetWaterTableInProfile()
            AC70_struc(n)%ac70(t)%StartMode = GetStartMode()
            AC70_struc(n)%ac70(t)%NoMoreCrop = GetNoMoreCrop()
            AC70_struc(n)%ac70(t)%CGCadjustmentAfterCutting = GetCGCadjustmentAfterCutting()
          
            ! calculate Aquacrop vegetation parameter for WCM and write as LAI
            ! as expected in WCM call (?)
            ! shift all previous Tmin by 1 index
            do itemp = AC70_struc(n)%Tmin_windowsize, 2, -1
                AC70_struc(n)%ac70(t)%Tmin_ac_antecedent(itemp) = AC70_struc(n)%ac70(t)%Tmin_ac_antecedent(itemp-1) ! degC
            end do
            ! add new Tmin on position 1
            AC70_struc(n)%ac70(t)%Tmin_ac_antecedent(1) = AC70_struc(n)%ac70(t)%Tmin_ac ! degC
            Tmin_movmean = 0.0
            countertemp = 0
            do itemp = 1, AC70_struc(n)%Tmin_windowsize
                Tmin_movmean = Tmin_movmean + AC70_struc(n)%ac70(t)%Tmin_ac_antecedent(itemp)
                countertemp = countertemp + 1
            end do 
            Tmin_movmean = Tmin_movmean / countertemp
            ! before 20 July
            if ((LIS_rc%mo .eq. 7) .and. (LIS_rc%da .eq. 20)) then
                AC70_struc(n)%ac70(t)%Tmin_ac_Julyref = Tmin_movmean
            end if
            if ((LIS_rc%mo .lt. 7) .or. ((LIS_rc%mo .eq. 7) .and. (LIS_rc%da .lt. 20))) then
                Tmin_mplr = 1.0
            else
                Tmin_mplr = Tmin_movmean / AC70_struc(n)%ac70(t)%Tmin_ac_Julyref
                Tmin_mplr = AMAX1(AMIN1(1.0,Tmin_mplr),0.0)
            end if
            if ((AC70_struc(n)%ac70(t)%SumWaBal%Biomass * Tmin_mplr) .gt. 0.0) then
                AC70_struc(n)%ac70(t)%WCMV1V2 = (AC70_struc(n)%ac70(t)%SumWaBal%Biomass * Tmin_mplr)**0.5
            else
                AC70_struc(n)%ac70(t)%WCMV1V2 = 0.0
            end if
            AC70_struc(n)%ac70(t)%AC70FC = AC70_struc(n)%ac70(t)%SoilLayer(1)%fc

    if ((LIS_rc%mo .eq. 12) .AND. (LIS_rc%da .eq. 31)) then
        AC70_struc(n)%ac70(t)%InitializeRun = 1
        !call FinalizeRun1(AC70_struc(n)%ac70(t)%irun, GetTheProjectFile(), AC70_struc(n)%ac70(t)%TheProjectType)
        !call FinalizeRun2(AC70_struc(n)%ac70(t)%irun, AC70_struc(n)%ac70(t)%TheProjectType)
        AC70_struc(n)%ac70(t)%irun = AC70_struc(n)%ac70(t)%irun + 1
    end if
            !!! MB_AC70


            ! MB: AC70
            ![ 1] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture, ice + liquid 
            do i=1, AC70_struc(n)%ac70(t)%NrCompartments
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC70_struc(n)%ac70(t)%smc(i),  &
                                                  vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC70_struc(n)%ac70(t)%smc(i), &
                                                  vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            ![ 4] output variable: biomass (unit=t/ha). ***  leaf area index 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70BIOMASS, value = real(AC70_struc(n)%ac70(t)%SumWaBal%Biomass,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70BIOMASS, value = real(AC70_struc(n)%ac70(t)%SumWaBal%Biomass,kind=sp), &
            !                                  vlevel=1, unit="t h-1", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CCiPrev, value = real(AC70_struc(n)%ac70(t)%CCiPrev,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_WCMV1V2, value = real(AC70_struc(n)%ac70(t)%WCMV1V2,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70FC, value = real(AC70_struc(n)%ac70(t)%AC70FC,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70Irrigation, value = real(AC70_struc(n)%ac70(t)%Irrigation,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70StressTot_Temp, value = real(AC70_struc(n)%ac70(t)%StressTot%Temp,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70StressTot_Sto, value = real(AC70_struc(n)%ac70(t)%StressTot%Sto,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70RootZoneWC_Actual, value = real(AC70_struc(n)%ac70(t)%RootZoneWC_Actual,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70RootZoneWC_WP, value = real(AC70_struc(n)%ac70(t)%RootZoneWC_WP,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70RootZoneWC_FC, value = real(AC70_struc(n)%ac70(t)%RootZoneWC_FC,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70Tact, value = real(AC70_struc(n)%ac70(t)%Tact,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70Eact, value = real(AC70_struc(n)%ac70(t)%Eact,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70RootingDepth, value = real(AC70_struc(n)%ac70(t)%RootingDepth,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC70CCiActual, value = real(AC70_struc(n)%ac70(t)%CCiActual,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PREC_ac, value = real(AC70_struc(n)%ac70(t)%PREC_ac,kind=sp), &
            !                                  vlevel=1, unit="mm/d", direction="-", surface_type = LIS_rc%lsm_index)
            
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PREC_ac, value=tmp_PREC_ac, &
            !    vlevel=1,unit="kg m-2 s-1",direction="DN",surface_type=LIS_rc%lsm_index)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_PREC_ac, value=tmp_PREC_ac*dt, &
            !    vlevel=1,unit="kg m-2", direction="DN",surface_type=LIS_rc%lsm_index)
            ! reset forcing variables to 0 for accumulation 

            AC70_struc(n)%ac70(t)%PREC_ac = 0.0
            AC70_struc(n)%ac70(t)%TMIN_ac = 0.0
            AC70_struc(n)%ac70(t)%TMAX_ac = 0.0
            AC70_struc(n)%ac70(t)%ETo_ac = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        AC70_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 

end subroutine Ac70_main
