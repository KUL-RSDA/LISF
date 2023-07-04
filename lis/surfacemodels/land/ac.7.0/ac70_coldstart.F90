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
! !ROUTINE: Ac70_coldstart
! \label{Ac70_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   9/4/14: Shugong Wang; initial implementation for LIS 7 and Ac70
!
! !INTERFACE:
subroutine Ac70_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use Ac70_lsmMod
 
!!! MB: AC70
    use ac_global, only: GetSimulParam_ThicknessTopSWC,&
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
                         GetTotalSaltContent,&
                         GetTotalWaterContent,&
                         Geteffectiverain,&
                         GetSumWaBal,&
                         GetRootZoneSalt,&
                         GetSimulation,&
                         GetCompartment,&
                         GetCompartment_theta,&
                         GetSoilLayer,&
                         GetIrrigation,&
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
                         GetIrriBeforeSeason,&
                         GetIrriAfterSeason


    use ac_run, only: GetIrriInterval,&
                      GetIrriInfoRecord1,&
                      GetIrriInfoRecord2,&

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
                    GetDayNri
              
              
    use ac_kinds, only: dp
!
! !DESCRIPTION:
!
!  This routine initializes the Ac70 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    
    ! added by Shugong Wang
    integer :: isnow
    real, allocatable, dimension(:) :: zsnso 
    real, allocatable, dimension(:) :: tsno
    real, allocatable, dimension(:) :: snice
    real, allocatable, dimension(:) :: snliq 
    real, allocatable, dimension(:) :: zsoil 
    ! end add

    !EMK...Temporary arrays.
    real :: tmp_swe(1,1),tmp_tgxy(1,1),tmp_snodep(1,1)
    integer :: tmp_isnowxy(1,1)

    do n=1, LIS_rc%nnest
       isnow = -AC70_struc(n)%nsnow
        ! added by shugong
        allocate(zsnso(-AC70_struc(n)%nsnow+1:AC70_struc(n)%nsoil))
        allocate(tsno(-AC70_struc(n)%nsnow+1:0))
        allocate(snice(-AC70_struc(n)%nsnow+1:0))
        allocate(snliq(-AC70_struc(n)%nsnow+1:0))
        allocate(zsoil(AC70_struc(n)%nsoil))
        zsoil(1) = -AC70_struc(n)%sldpth(1)
        do l=2, AC70_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - AC70_struc(n)%sldpth(l) 
        enddo
        ! end add 

        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Ac70_coldstart -- cold-starting Ac70"
            do t=1, LIS_rc%npatch(n,mtype)
                AC70_struc(n)%ac70(t)%albold = AC70_struc(n)%init_albold
                AC70_struc(n)%ac70(t)%sneqvo = AC70_struc(n)%init_sneqvo
                ! only soil temperature is intialized, snow temperature is calculated by snow_init 
                do l=1, AC70_struc(n)%nsoil
                    AC70_struc(n)%ac70(t)%sstc(AC70_struc(n)%nsnow+l) = AC70_struc(n)%init_stc(l)
                enddo
                do l=1, AC70_struc(n)%nsoil
                    AC70_struc(n)%ac70(t)%sh2o(l) = AC70_struc(n)%init_sh2o(l)
                enddo
                AC70_struc(n)%ac70(t)%tah = AC70_struc(n)%init_tah
                AC70_struc(n)%ac70(t)%eah = AC70_struc(n)%init_eah
                AC70_struc(n)%ac70(t)%fwet = AC70_struc(n)%init_fwet
                AC70_struc(n)%ac70(t)%canliq = AC70_struc(n)%init_canliq
                AC70_struc(n)%ac70(t)%canice = AC70_struc(n)%init_canice
                AC70_struc(n)%ac70(t)%tv = AC70_struc(n)%init_tv
                AC70_struc(n)%ac70(t)%tg = AC70_struc(n)%init_tg
                AC70_struc(n)%ac70(t)%qsnow = AC70_struc(n)%init_qsnow
                !do l=1, AC70_struc(n)%nsoil + AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%zss(l) = AC70_struc(n)%init_zss(l)
                !enddo
                AC70_struc(n)%ac70(t)%snowh = AC70_struc(n)%init_snowh
                AC70_struc(n)%ac70(t)%sneqv = AC70_struc(n)%init_sneqv
                !do l=1, AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%snowice(l) = AC70_struc(n)%init_snowice(l)
                !enddo
                !do l=1, AC70_struc(n)%nsnow
                !    AC70_struc(n)%ac70(t)%snowliq(l) = AC70_struc(n)%init_snowliq(l)
                !enddo
                AC70_struc(n)%ac70(t)%zwt = AC70_struc(n)%init_zwt
                AC70_struc(n)%ac70(t)%wa = AC70_struc(n)%init_wa
                AC70_struc(n)%ac70(t)%wt = AC70_struc(n)%init_wt
                AC70_struc(n)%ac70(t)%wslake = AC70_struc(n)%init_wslake
                AC70_struc(n)%ac70(t)%lfmass = AC70_struc(n)%init_lfmass
                AC70_struc(n)%ac70(t)%rtmass = AC70_struc(n)%init_rtmass
                AC70_struc(n)%ac70(t)%stmass = AC70_struc(n)%init_stmass
                AC70_struc(n)%ac70(t)%wood = AC70_struc(n)%init_wood
                AC70_struc(n)%ac70(t)%stblcp = AC70_struc(n)%init_stblcp
                AC70_struc(n)%ac70(t)%fastcp = AC70_struc(n)%init_fastcp
                AC70_struc(n)%ac70(t)%lai = AC70_struc(n)%init_lai
                AC70_struc(n)%ac70(t)%sai = AC70_struc(n)%init_sai
                AC70_struc(n)%ac70(t)%cm = AC70_struc(n)%init_cm
                AC70_struc(n)%ac70(t)%ch = AC70_struc(n)%init_ch
                AC70_struc(n)%ac70(t)%tauss = AC70_struc(n)%init_tauss
                AC70_struc(n)%ac70(t)%smcwtd = AC70_struc(n)%init_smcwtd
                AC70_struc(n)%ac70(t)%deeprech = AC70_struc(n)%init_deeprech
                AC70_struc(n)%ac70(t)%rech = AC70_struc(n)%init_rech
                AC70_struc(n)%ac70(t)%zlvl = AC70_struc(n)%init_zlvl 
                !!! MB: AC70
                !AC70_struc(n)%ac70(t)%daynri = AC70_struc(n)%init_daynri
                !AC70_struc(n)%ac70(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
                !AC70_struc(n)%ac70(t)%RootZoneWC_FC = GetRootZoneWC_FC()
                !AC70_struc(n)%ac70(t)%RootZoneWC_WP = GetRootZoneWC_WP()
                !AC70_struc(n)%ac70(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
                !AC70_struc(n)%ac70(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
                !AC70_struc(n)%ac70(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
                !AC70_struc(n)%ac70(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
                !AC70_struc(n)%ac70(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
                !AC70_struc(n)%ac70(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
                !AC70_struc(n)%ac70(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
                !AC70_struc(n)%ac70(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
                !AC70_struc(n)%ac70(t)%Compartment = GetCompartment()
                !AC70_struc(n)%ac70(t)%soillayer = GetSoilLayer()
                !AC70_struc(n)%ac70(t)%TotalSaltContent = GetTotalSaltContent()
                !AC70_struc(n)%ac70(t)%TotalWaterContent = GetTotalWaterContent()
                !AC70_struc(n)%ac70(t)%effectiverain = Geteffectiverain()
                !AC70_struc(n)%ac70(t)%SumWaBal = GetSumWaBal()
                !AC70_struc(n)%ac70(t)%RootZoneSalt = GetRootZoneSalt()
                !AC70_struc(n)%ac70(t)%Simulation = GetSimulation()
                !AC70_struc(n)%ac70(t)%IrriInterval = GetIrriInterval()
                !AC70_struc(n)%ac70(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
                !AC70_struc(n)%ac70(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
                !AC70_struc(n)%ac70(t)%Irrigation = GetIrrigation()
                !AC70_struc(n)%ac70(t)%IrriBeforeSeason = GetIrriBeforeSeason()
                !AC70_struc(n)%ac70(t)%IrriAfterSeason = GetIrriAfterSeason()
                !AC70_struc(n)%ac70(t)%IrriECw = GetIrriECw()
                !AC70_struc(n)%ac70(t)%Management = GetManagement()
                !AC70_struc(n)%ac70(t)%PerennialPeriod = GetPerennialPeriod()
                !AC70_struc(n)%ac70(t)%simulparam = GetSimulParam()
                !AC70_struc(n)%ac70(t)%Cuttings = GetManagement_Cuttings()
                !AC70_struc(n)%ac70(t)%onset = GetOnset()
                !AC70_struc(n)%ac70(t)%endseason = GetEndSeason()
                !AC70_struc(n)%ac70(t)%crop = GetCrop()
                !AC70_struc(n)%ac70(t)%Soil = GetSoil()
                !AC70_struc(n)%ac70(t)%TemperatureRecord = GetTemperatureRecord()
                !AC70_struc(n)%ac70(t)%ClimRecord = GetClimRecord()
                !AC70_struc(n)%ac70(t)%RainRecord = GetRainRecord()
                !AC70_struc(n)%ac70(t)%EToRecord = GetEToRecord()
                !AC70_struc(n)%ac70(t)%daynri = GetDayNri()
                AC70_struc(n)%ac70(t)%irun = 1
                !AC70_struc(n)%daynrinextclimaterecord = GetDayNri() + 1
                AC70_struc(n)%daynrinextclimaterecord = 1

                !AC70_struc(n)%ac70(t)%GenerateTimeMode = GetGenerateTimeMode()
                !AC70_struc(n)%ac70(t)%GenerateDepthMode = GetGenerateDepthMode()
                !AC70_struc(n)%ac70(t)%IrriMode = GetIrriMode()
                !AC70_struc(n)%ac70(t)%IrriMethod = GetIrriMethod()
                !AC70_struc(n)%ac70(t)%DaySubmerged = GetDaySubmerged()
                !AC70_struc(n)%ac70(t)%MaxPlotNew = GetMaxPlotNew()
                !AC70_struc(n)%ac70(t)%NrCompartments = GetNrCompartments()
                !AC70_struc(n)%ac70(t)%IrriFirstDayNr = GetIrriFirstDayNr()
                !AC70_struc(n)%ac70(t)%ZiAqua = GetZiAqua()
                !AC70_struc(n)%ac70(t)%IniPercTAW = GetIniPercTAW()
                !AC70_struc(n)%ac70(t)%MaxPlotTr = GetMaxPlotTr()
                !AC70_struc(n)%ac70(t)%OutputAggregate = GetOutputAggregate()

                !AC70_struc(n)%ac70(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
                !AC70_struc(n)%ac70(t)%PreDay = GetPreDay()
                !AC70_struc(n)%ac70(t)%OutDaily = GetOutDaily()
                !AC70_struc(n)%ac70(t)%Out1Wabal = GetOut1Wabal()
                !AC70_struc(n)%ac70(t)%Out2Crop = GetOut2Crop()
                !AC70_struc(n)%ac70(t)%Out3Prof = GetOut3Prof()
                !AC70_struc(n)%ac70(t)%Out4Salt = GetOut4Salt()
                !AC70_struc(n)%ac70(t)%Out5CompWC = GetOut5CompWC()
                !AC70_struc(n)%ac70(t)%Out6CompEC = GetOut6CompEC()
                !AC70_struc(n)%ac70(t)%Out7Clim = GetOut7Clim()
                !AC70_struc(n)%ac70(t)%Part1Mult = GetPart1Mult()
                !AC70_struc(n)%ac70(t)%Part2Eval = GetPart2Eval()

                !
                !AC70_struc(n)%ac70(t)%CCiActual = GetCCiActual()
                !AC70_struc(n)%ac70(t)%CCiprev = GetCCiprev()
                !AC70_struc(n)%ac70(t)%CCiTopEarlySen = GetCCiTopEarlySen()
                !AC70_struc(n)%ac70(t)%CRsalt = GetCRsalt () ! gram/m2
                !AC70_struc(n)%ac70(t)%CRwater = GetCRwater() ! mm/day
                !AC70_struc(n)%ac70(t)%ECdrain = GetECdrain() ! EC drain water dS/m
                !AC70_struc(n)%ac70(t)%ECiAqua = GetECiAqua() ! EC of the groundwater table in dS/m
                !AC70_struc(n)%ac70(t)%ECstorage = GetECstorage() !EC surface storage dS/m
                !AC70_struc(n)%ac70(t)%Eact = GetEact() ! mm/day
                !AC70_struc(n)%ac70(t)%Epot = GetEpot() ! mm/day
                !AC70_struc(n)%ac70(t)%Drain = GetDrain()  ! mm/day
                !AC70_struc(n)%ac70(t)%Infiltrated = GetInfiltrated() ! mm/day
                !AC70_struc(n)%ac70(t)%RootingDepth = GetRootingDepth()
                !AC70_struc(n)%ac70(t)%Runoff = GetRunoff()  ! mm/day
                !AC70_struc(n)%ac70(t)%SaltInfiltr = GetSaltInfiltr() ! salt infiltrated in soil profile Mg/ha
                !AC70_struc(n)%ac70(t)%Surf0 = GetSurf0()  ! surface water [mm] begin day
                !AC70_struc(n)%ac70(t)%SurfaceStorage = GetSurfaceStorage() !mm/day
                !AC70_struc(n)%ac70(t)%Tact = GetTact() ! mm/day
                !AC70_struc(n)%ac70(t)%Tpot = GetTpot() ! mm/day
                !AC70_struc(n)%ac70(t)%TactWeedInfested = GetTactWeedInfested() !mm/day

                !if (trim(LIS_rc%metforc(1)) == 'MERRA2_AC') then
                !   AC70_struc(n)%ac70(t)%PREC_ac = 0.0  ! mm/day
                !   AC70_struc(n)%ac70(t)%Tmin_ac = 0.0 ! degC
                !   AC70_struc(n)%ac70(t)%Tmax_ac = 0.0 ! degC
                !   AC70_struc(n)%ac70(t)%ETo_ac = 0.0 ! mm/day
                !else
                !   AC70_struc(n)%ac70(t)%PREC_ac = GetRain()  ! mm/day
                !   AC70_struc(n)%ac70(t)%Tmin_ac =GetTmin() ! degC
                !   AC70_struc(n)%ac70(t)%Tmax_ac = GetTmax() ! degC
                !   AC70_struc(n)%ac70(t)%ETo_ac = GetETo() ! mm/day
                !end if

                AC70_struc(n)%ac70(t)%InitializeRun = 1 ! gets 1 at end of year 

                !AC70_struc(n)%ac70(t)%GwTable = GetGwTable()
                !AC70_struc(n)%ac70(t)%PlotVarCrop = GetPlotVarCrop()
                !AC70_struc(n)%ac70(t)%StressTot = GetStressTot()
                !AC70_struc(n)%ac70(t)%CutInfoRecord1 = GetCutInfoRecord1()
                !AC70_struc(n)%ac70(t)%CutInfoRecord2 = GetCutInfoRecord2()
                !AC70_struc(n)%ac70(t)%Transfer = GetTransfer()
                !AC70_struc(n)%ac70(t)%PreviousSum = GetPreviousSum()
                !AC70_struc(n)%ac70(t)%Tadj = GetTadj()
                !AC70_struc(n)%ac70(t)%GDDTadj = GetGDDTadj()
                !AC70_struc(n)%ac70(t)%DayLastCut = GetDayLastCut()
                !AC70_struc(n)%ac70(t)%NrCut = GetNrCut()
                !AC70_struc(n)%ac70(t)%SumInterval = GetSumInterval()
                !AC70_struc(n)%ac70(t)%PreviousStressLevel = GetPreviousStressLevel()
                !AC70_struc(n)%ac70(t)%StressSFadjNEW = GetStressSFadjNEW()
                !AC70_struc(n)%ac70(t)%Bin = GetBin()
                !AC70_struc(n)%ac70(t)%Bout = GetBout()
                !AC70_struc(n)%ac70(t)%GDDayi = GetGDDayi()
                !AC70_struc(n)%ac70(t)%CO2i = GetCO2i()
                !AC70_struc(n)%ac70(t)%FracBiomassPotSF = GetFracBiomassPotSF()
                !AC70_struc(n)%ac70(t)%SumETo = GetSumETo()
                !AC70_struc(n)%ac70(t)%SumGDD = GetSumGDD()
                !AC70_struc(n)%ac70(t)%Ziprev = GetZiprev()
                !AC70_struc(n)%ac70(t)%SumGDDPrev = GetSumGDDPrev()
                !AC70_struc(n)%ac70(t)%CCxWitheredTpot = GetCCxWitheredTpot()
                !AC70_struc(n)%ac70(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
                !AC70_struc(n)%ac70(t)%Coeffb0 = GetCoeffb0()
                !AC70_struc(n)%ac70(t)%Coeffb1 = GetCoeffb1()
                !AC70_struc(n)%ac70(t)%Coeffb2 = GetCoeffb2()
                !AC70_struc(n)%ac70(t)%Coeffb0Salt = GetCoeffb0Salt()
                !AC70_struc(n)%ac70(t)%Coeffb1Salt = GetCoeffb1Salt()
                !AC70_struc(n)%ac70(t)%Coeffb2Salt = GetCoeffb2Salt()
                !AC70_struc(n)%ac70(t)%StressLeaf = GetStressLeaf()
                !AC70_struc(n)%ac70(t)%StressSenescence = GetStressSenescence ()
                !AC70_struc(n)%ac70(t)%DayFraction = GetDayFraction()
                !AC70_struc(n)%ac70(t)%GDDayFraction = GetGDDayFraction()
                !AC70_struc(n)%ac70(t)%CGCref = GetCGCref()
                !AC70_struc(n)%ac70(t)%GDDCGCref = GetGDDCGCref ()
                !AC70_struc(n)%ac70(t)%TimeSenescence = GetTimeSenescence ()
                !AC70_struc(n)%ac70(t)%SumKcTop = GetSumKcTop()
                !AC70_struc(n)%ac70(t)%SumKcTopStress = GetSumKcTopStress()
                !AC70_struc(n)%ac70(t)%SumKci = GetSumKci()
                !AC70_struc(n)%ac70(t)%CCoTotal = GetCCoTotal()
                !AC70_struc(n)%ac70(t)%CCxTotal = GetCCxTotal()
                !AC70_struc(n)%ac70(t)%CDCTotal = GetCDCTotal()
                !AC70_struc(n)%ac70(t)%GDDCDCTotal = GetGDDCDCTotal()
                !AC70_struc(n)%ac70(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
                !AC70_struc(n)%ac70(t)%WeedRCi = GetWeedRCi()
                !AC70_struc(n)%ac70(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
                !AC70_struc(n)%ac70(t)%fWeedNoS = GetfWeedNoS()
                !AC70_struc(n)%ac70(t)%Zeval = GetZeval()
                !AC70_struc(n)%ac70(t)%BprevSum = GetBprevSum()
                !AC70_struc(n)%ac70(t)%YprevSum = GetYprevSum()
                !AC70_struc(n)%ac70(t)%SumGDDcuts = GetSumGDDcuts()
                !AC70_struc(n)%ac70(t)%HItimesBEF = GetHItimesBEF()
                !AC70_struc(n)%ac70(t)%ScorAT1 = GetScorAT1()
                !AC70_struc(n)%ac70(t)%ScorAT2 = GetScorAT2()
                !AC70_struc(n)%ac70(t)%HItimesAT1 = GetHItimesAT1()
                !AC70_struc(n)%ac70(t)%HItimesAT2 = GetHItimesAT2()
                !AC70_struc(n)%ac70(t)%HItimesAT = GetHItimesAT()
                !AC70_struc(n)%ac70(t)%alfaHI = GetalfaHI()
                !AC70_struc(n)%ac70(t)%alfaHIAdj = GetalfaHIAdj()
                !AC70_struc(n)%ac70(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
                !AC70_struc(n)%ac70(t)%DayNr1Eval = GetDayNr1Eval()
                !AC70_struc(n)%ac70(t)%DayNrEval = GetDayNrEval()
                !AC70_struc(n)%ac70(t)%LineNrEval = GetLineNrEval()
                !AC70_struc(n)%ac70(t)%PreviousSumETo = GetPreviousSumETo()
                !AC70_struc(n)%ac70(t)%PreviousSumGDD = GetPreviousSumGDD()
                !AC70_struc(n)%ac70(t)%PreviousBmob = GetPreviousBmob()
                !AC70_struc(n)%ac70(t)%PreviousBsto = GetPreviousBsto()
                !AC70_struc(n)%ac70(t)%StageCode = GetStageCode()
                !AC70_struc(n)%ac70(t)%PreviousDayNr = GetPreviousDayNr()
                !AC70_struc(n)%ac70(t)%NoYear = GetNoYear()
                !AC70_struc(n)%ac70(t)%WaterTableInProfile = GetWaterTableInProfile()
                !AC70_struc(n)%ac70(t)%StartMode = GetStartMode()
                !AC70_struc(n)%ac70(t)%NoMoreCrop = GetNoMoreCrop()
                !AC70_struc(n)%ac70(t)%CGCadjustmentAfterCutting = GetCGCadjustmentAfterCutting()

                do l=1, AC70_struc(n)%ac70(t)%NrCompartments
                    AC70_struc(n)%ac70(t)%smc(l) = GetCompartment_theta(l)
                enddo

                do l=1, 40
                    AC70_struc(n)%ac70(t)%Tmin_ac_antecedent(l) = 0.0
                enddo
                ! added by shugong 
                zsnso = 0.0 
                !EMK...snow_init_70 is expecting several arrays which
                !are being passed as scalars.  Although no memory corruption
                !occurs here because of the declared array dimensions (all 1),
                !this is still technically a syntax error.  So, we will
                !copy the required fields to temporary arrays with the
                !correct declarations and pass those instead.
                !call snow_init_70(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                !                  AC70_struc(n)%nsnow,          & !input 
                !                  AC70_struc(n)%nsoil,          & !input 
                !                  zsoil,                            & !input
                !                  AC70_struc(n)%init_sneqv,     & !input
                !                  AC70_struc(n)%init_tg,        & !input
                !                  AC70_struc(n)%init_snowh,     & !input
                !                  zsnso, tsno, snice, snliq, isnow) ! output 
                tmp_swe(1,1) = AC70_struc(n)%init_sneqv
                tmp_tgxy(1,1) = AC70_struc(n)%init_tg
                tmp_snodep(1,1) = AC70_struc(n)%init_snowh
                call snow_init_70(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                                  AC70_struc(n)%nsnow,          & !input 
                                  AC70_struc(n)%nsoil,          & !input 
                                  zsoil,                            & !input
                                  tmp_swe,         & !input
                                  tmp_tgxy,        & !input
                                  tmp_snodep,      & !input
                                  zsnso, tsno, snice, snliq, &
                                  tmp_isnowxy) ! output
                isnow = tmp_isnowxy(1,1)
                AC70_struc(n)%ac70(t)%snowice(1:AC70_struc(n)%nsnow) = snice(-AC70_struc(n)%nsnow+1:0)
                AC70_struc(n)%ac70(t)%snowliq(1:AC70_struc(n)%nsnow) = snliq(-AC70_struc(n)%nsnow+1:0)
                AC70_struc(n)%ac70(t)%zss(1:AC70_struc(n)%nsnow+AC70_struc(n)%nsoil) = zsnso(-AC70_struc(n)%nsnow+1:AC70_struc(n)%nsoil) 
                AC70_struc(n)%ac70(t)%sstc(AC70_struc(n)%nsnow+isnow+1:AC70_struc(n)%nsnow) = tsno(isnow+1:0) 
                AC70_struc(n)%ac70(t)%isnow = isnow                

                
                ! end add
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: Ac70_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
        deallocate(zsnso)
        deallocate(tsno)
        deallocate(snice)
        deallocate(snliq)
        deallocate(zsoil) 
    enddo
end subroutine Ac70_coldstart
