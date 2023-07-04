!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: Ac70_set_pedecvars
!  \label{Ac70_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!


! !INTERFACE:
subroutine Ac70_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use Ac70_lsmMod, only : Ac70_struc
  use Ac70_peMod,  only : Ac70_pe_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
! 
!EOP
  integer                :: n
  real, pointer          :: vdata(:)
  character*100          :: vname
  integer, pointer       :: mod_flag_Ac70(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_Ac70(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_Ac70 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,Ac70_pe_struc(n)%nparams
     if(Ac70_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(Ac70_pe_struc(n)%param_name(i))
        call Ac70_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call Ac70_checkBounds(n,DEC_State,vname, vdata, mod_flag_Ac70)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call Ac70_checkConstraints(n,DEC_State, mod_flag_Ac70)

  !set variables given modflag; if flag set will leave values alone
  call Ac70_setVars(n,DEC_State,mod_flag_Ac70)

  !send mod flag to ESMF state (feasibility flag)
  call Ac70_setModFlag(n,DEC_State,Feas_State,mod_flag_Ac70)
end subroutine Ac70_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine Ac70_getvardata(n,DEC_State,vname, vdata, statusStateGet)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vdata(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
! 
!EOP
  real, pointer          :: vardata(:)
  type(ESMF_Field)       :: varField
  integer                :: statusStateGet, statusFieldGet,i
  
  call ESMF_StateGet(DEC_State,vname,varField,rc=statusStateGet)
!  call LIS_verify(status)
  
  if(statusStateGet.eq.0) then
     call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
          rc=statusFieldGet)
     call LIS_verify(statusFieldGet)
     vdata=vardata
  endif
  
end subroutine Ac70_getvardata

subroutine Ac70_checkBounds(n,DEC_State,vname, vardata, mod_flag_Ac70)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer       :: mod_flag_Ac70(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  integer                :: status
  integer                :: t

  call ESMF_StateGet(DEC_State,vname,varField,rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata(t).lt.vardata_min) then 
        mod_flag_Ac70(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_Ac70(t) = 1
     endif
  enddo
end subroutine Ac70_checkBounds

subroutine Ac70_checkConstraints(n,DEC_State,mod_flag_Ac70)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use Ac70_lsmMod, only : Ac70_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_Ac70(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  character*100          :: vname
  integer                :: t
  integer       :: status1, status2
  real, allocatable :: vardata1(:)
  real, allocatable :: vardata2(:)

  allocate(vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  !SMCMAX > SMCDRY
  vname='SMCMAX'
  call Ac70_getvardata(n,DEC_State,vname,vardata1, status1)
!  vname='SMCDRY'
!  call Ac70_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=Ac70_struc(n)%ac70(:)%smcmax
!  if(status2.ne.0) vardata2=Ac70_struc(n)%ac70(:)%smcdry
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(vardata1(t).le.vardata2(t)) then
     if(vardata1(t).le.Ac70_struc(n)%ac70(t)%smcdry) then
        mod_flag_Ac70(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call Ac70_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call Ac70_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=Ac70_struc(n)%ac70(:)%smcref
  if(status2.ne.0) vardata2=Ac70_struc(n)%ac70(:)%smcwlt
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_Ac70(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call Ac70_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call Ac70_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=Ac70_struc(n)%ac70(:)%smcmax
  if(status2.ne.0) vardata2=Ac70_struc(n)%ac70(:)%smcref
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_Ac70(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call Ac70_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call Ac70_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=Ac70_struc(n)%ac70(:)%HVT
  if(status2.ne.0) vardata2=Ac70_struc(n)%ac70(:)%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then ! SY: Note .lt. instead of .le., following some entries with HVT=HVB in MPTABLE_UMD.TBL
        mod_flag_Ac70(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call Ac70_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call Ac70_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=Ac70_struc(n)%ac70(:)%HVT
  if(status2.ne.0) vardata2=Ac70_struc(n)%ac70(:)%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then 
        mod_flag_Ac70(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)

end subroutine Ac70_checkConstraints

subroutine Ac70_setVars(n,DEC_State,mod_flag_Ac70)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use Ac70_lsmMod, only : Ac70_struc
  use Ac70_peMod,  only : Ac70_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_Ac70(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Ac70 model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,Ac70_pe_struc(n)%nparams
     if(Ac70_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(Ac70_pe_struc(n)%param_name(i))
        call Ac70_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_Ac70(t).eq.0) then 
              if(vname.eq."TOPT") &
                   Ac70_struc(n)%ac70(t)%topt = vardata(t) 
              if(vname.eq."RGL") & 
                   Ac70_struc(n)%ac70(t)%rgl = vardata(t) 
              if(vname.eq."RSMAX") &
                   Ac70_struc(n)%ac70(t)%rsmax = vardata(t) 
              if(vname.eq."RSMIN") & 
                   Ac70_struc(n)%ac70(t)%rsmin = vardata(t) 
              if(vname.eq."HS") & 
                   Ac70_struc(n)%ac70(t)%hs = vardata(t) 
              if(vname.eq."NROOT") &
                   Ac70_struc(n)%ac70(t)%nroot = vardata(t) 
              if(vname.eq."CSOIL") &
                   Ac70_struc(n)%ac70(t)%csoil = vardata(t) 
              if(vname.eq."BEXP") & 
                   Ac70_struc(n)%ac70(t)%bexp = vardata(t) 
              if(vname.eq."DKSAT") & 
                   Ac70_struc(n)%ac70(t)%dksat = vardata(t) 
              if(vname.eq."DWSAT") & 
                   Ac70_struc(n)%ac70(t)%dwsat = vardata(t) 
              if(vname.eq."PSISAT") & 
                   Ac70_struc(n)%ac70(t)%psisat = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   Ac70_struc(n)%ac70(t)%quartz = vardata(t) 
              if(vname.eq."SMCMAX") &
                   Ac70_struc(n)%ac70(t)%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   Ac70_struc(n)%ac70(t)%smcref = vardata(t) 
              if(vname.eq."SMCWLT") &
                   Ac70_struc(n)%ac70(t)%smcwlt = vardata(t) 
              if(vname.eq."CZIL") &
                   Ac70_struc(n)%ac70(t)%czil = vardata(t) 
              if(vname.eq."FRZK") &
                   Ac70_struc(n)%ac70(t)%frzk = vardata(t) 
              if(vname.eq."REFDK") &
                   Ac70_struc(n)%ac70(t)%refdk = vardata(t) 
              if(vname.eq."REFKDT") &
                   Ac70_struc(n)%ac70(t)%refkdt = vardata(t) 
              if(vname.eq."SLOPE") &
                   Ac70_struc(n)%ac70(t)%slope = vardata(t) 
              if(vname.eq."CH2OP") &
                   Ac70_struc(n)%ac70(t)%CH2OP = vardata(t)
              if(vname.eq."DLEAF") &
                   Ac70_struc(n)%ac70(t)%DLEAF = vardata(t)
              if(vname.eq."Z0MVT") &
                   Ac70_struc(n)%ac70(t)%Z0MVT = vardata(t)
              if(vname.eq."HVT") &
                   Ac70_struc(n)%ac70(t)%HVT = vardata(t)
              if(vname.eq."HVB") &
                   Ac70_struc(n)%ac70(t)%HVB = vardata(t)
              if(vname.eq."RC") &
                   Ac70_struc(n)%ac70(t)%RC = vardata(t)
              if(vname.eq."RHOL1") &
                   Ac70_struc(n)%ac70(t)%RHOL1 = vardata(t)
              if(vname.eq."RHOL2") &
                   Ac70_struc(n)%ac70(t)%RHOL2 = vardata(t)
              if(vname.eq."RHOS1") &
                   Ac70_struc(n)%ac70(t)%RHOS1 = vardata(t)
              if(vname.eq."RHOS2") &
                   Ac70_struc(n)%ac70(t)%RHOS2 = vardata(t)
              if(vname.eq."TAUL1") &
                   Ac70_struc(n)%ac70(t)%TAUL1 = vardata(t)
              if(vname.eq."TAUL2") &
                   Ac70_struc(n)%ac70(t)%TAUL2 = vardata(t)
              if(vname.eq."TAUS1") &
                   Ac70_struc(n)%ac70(t)%TAUS1 = vardata(t)
              if(vname.eq."TAUS2") &
                   Ac70_struc(n)%ac70(t)%TAUS2 = vardata(t)
              if(vname.eq."XL") &
                   Ac70_struc(n)%ac70(t)%XL = vardata(t)
              if(vname.eq."CWPVT") &
                   Ac70_struc(n)%ac70(t)%CWPVT = vardata(t)
              if(vname.eq."C3PSN") &
                   Ac70_struc(n)%ac70(t)%C3PSN = vardata(t)
              if(vname.eq."KC25") &
                   Ac70_struc(n)%ac70(t)%KC25 = vardata(t)
              if(vname.eq."AKC") &
                   Ac70_struc(n)%ac70(t)%AKC = vardata(t)
              if(vname.eq."KO25") &
                   Ac70_struc(n)%ac70(t)%KO25 = vardata(t)
              if(vname.eq."AKO") &
                   Ac70_struc(n)%ac70(t)%AKO = vardata(t)
              if(vname.eq."AVCMX") &
                   Ac70_struc(n)%ac70(t)%AVCMX = vardata(t)
              if(vname.eq."AQE") &
                   Ac70_struc(n)%ac70(t)%AQE = vardata(t)
              if(vname.eq."LTOVRC") &
                   Ac70_struc(n)%ac70(t)%LTOVRC = vardata(t)
              if(vname.eq."DILEFC") &
                   Ac70_struc(n)%ac70(t)%DILEFC = vardata(t)
              if(vname.eq."DILEFW") &
                   Ac70_struc(n)%ac70(t)%DILEFW = vardata(t)
              if(vname.eq."RMF25") &
                   Ac70_struc(n)%ac70(t)%RMF25 = vardata(t)
              if(vname.eq."SLA") &
                   Ac70_struc(n)%ac70(t)%SLA = vardata(t)
              if(vname.eq."FRAGR") &
                   Ac70_struc(n)%ac70(t)%FRAGR = vardata(t)
              if(vname.eq."TMIN") &
                   Ac70_struc(n)%ac70(t)%TMIN = vardata(t)
              if(vname.eq."VCMX25") &
                   Ac70_struc(n)%ac70(t)%VCMX25 = vardata(t)
              if(vname.eq."TDLEF") &
                   Ac70_struc(n)%ac70(t)%TDLEF = vardata(t)
              if(vname.eq."BP") &
                   Ac70_struc(n)%ac70(t)%BP = vardata(t)
              if(vname.eq."MP") &
                   Ac70_struc(n)%ac70(t)%MP = vardata(t)
              if(vname.eq."QE25") &
                   Ac70_struc(n)%ac70(t)%QE25 = vardata(t)
              if(vname.eq."RMS25") &
                   Ac70_struc(n)%ac70(t)%RMS25 = vardata(t)
              if(vname.eq."RMR25") &
                   Ac70_struc(n)%ac70(t)%RMR25 = vardata(t)
              if(vname.eq."ARM") &
                   Ac70_struc(n)%ac70(t)%ARM = vardata(t)
              if(vname.eq."FOLNMX") &
                   Ac70_struc(n)%ac70(t)%FOLNMX = vardata(t)
              if(vname.eq."WDPOOL") &
                   Ac70_struc(n)%ac70(t)%WDPOOL = vardata(t)
              if(vname.eq."WRRAT") &
                   Ac70_struc(n)%ac70(t)%WRRAT = vardata(t)
              if(vname.eq."MRP") &
                   Ac70_struc(n)%ac70(t)%MRP = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine Ac70_setVars

subroutine Ac70_setModFlag(n,DEC_State,Feas_State,mod_flag_Ac70)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_Ac70(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine sets the feasibility flag
! 
!EOP
  type(ESMF_Field)       :: feasField
  integer                :: t
  integer                :: status
  integer, pointer       :: modflag(:)

  call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField,localDE=0,farrayPtr=modflag,rc=status)
  call LIS_verify(status)
!  write(LIS_logunit,*) 'Ac70_setModFlag 1 modflag:', modflag

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(mod_flag_Ac70(t).eq.1) then 
        modflag(t)=1
     endif
  enddo
!  write(LIS_logunit,*) 'Ac70_setModFlag 2 modflag:', modflag
end subroutine Ac70_setModFlag
