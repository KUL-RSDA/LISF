C MEMBER FROST2
C REPLACEMENT OF OLD FROST1
C
      SUBROUTINE FROST2_1(PX,TA,WE,AESC,SH,FRZPAR,SACPAR,FRZST,SACST,
     +            SACST_PRV,SMC,SH2O,DTFRZ,IDT,NSOIL,NUPL,NSAC,IVERS,
c     +            FRZDUP,FRZDBT,FROST)
     +            FRZDUP,FRZDBT,FROST,SMAX, SMCDRY, HRAPX,HRAPY)
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    PURPOSE:  TO CALCULATE SOIL TEMPERATURE & ICE CONTENT
C     WRITTEN BY VICTOR KOREN - HRL   JANUARY 2000
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  DEFAULT CONSTANTS
      PARAMETER (T0 = 273.16)
      PARAMETER (FRST_FACT = 5.0)

C  VARIABLE IVERS CONTROLL FROZEN GROUND OPTION:
C  IVERS = 0 FULLY NON-FROZEN GROUND
C  IVERS = 1 FULLY FROZEN GROUND (NEW VERSION)
C  IVERS = 4 WILL RUN SOIL TEMPERATURE BUT NOT FREEZ SOIL MOISTURE
      
      REAL FRZPAR(*),SACPAR(*),FRZST(*),SACST(*),SACST_PRV(*)
      REAL SMC(*), SH2O(*)
      REAL SAND(12),CLAY(12)
      
C--      REAL SMC   ( 10 )
C--      REAL SH2O  ( 10 )
      REAL RHSTS ( 10 )
      REAL STC   ( 10 )
      REAL STCOUT ( 10 )
      REAL ZSOIL ( 10 )
C--      REAL DSW   ( 10 )
C--      INTEGER IDSW (10)

cc      SAVE NSOIL,ZSOIL,TBOT,ZBOT
C--      SAVE NSOIL,ZSOIL,TBOT
C--      SAVE UZTWC0,UZFWC0,LZTWC0,LZFSC0,LZFPC0
C--      SAVE SMC,SH2O
C--      SAVE DTFRZ,ISTART,IDT,NLOOP
      
C--      INTEGER ISTART/0/,NLOOP/0/
C--      REAL SMAX,PSISAT,BRT,QUARTZ,STYPE
C--      REAL DT,PX,TA,WE,AESC
      REAL LZTWM,LZFSM,LZFPM,LZSK,LZPK,LZTWC,LZFSC,LZFPC
      REAL LZTWH,LZFSH,LZFPH,LZTWC0,LZFSC0,LZFPC0

CVK_02  NEW COMMON STATEMENT FOR DESIRED SOIL LAYERS
C--      INTEGER NINT/5/,NINTW/5/
C--      REAL TSINT(10),SWINT(10),SWHINT(10)
C--      REAL DSINT(10)/0.05,0.1,0.2,0.5,1.0,1.5,2.0,0.,0.,0./
C--      REAL DSINTW(10)/0.05,0.1,0.2,0.5,1.0,1.5,2.0,0.,0.,0./
C--      REAL TSTMP(10),DSMOD(10),SWTMP(10),SWHTMP(10)
C--      SAVE DSINT,DSINTW
C--      COMMON/TSLINT/TSINT,NINT,SWINT,SWHINT,NINTW      

      DATA SAND/0.92,0.82,0.58,0.17,0.09,0.43,0.58,0.10,0.32,0.52,
     +          0.06,0.22/
      DATA CLAY/0.03,0.06,0.10,0.13,0.05,0.18,0.27,0.34,0.34,0.42,
     +          0.47,0.58/      

C     COMMON BLOCKS
C--      COMMON/IONUM/IN,IPR,IPU
cc      COMMON/FSMPM1/UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,
cc     1              LZFSM,LZFPM,LZSK,LZPK,PFREE,SIDE,SAVED,PAREA
C--      COMMON/FSMCO1/UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,FGCO(6),RSUM(7),
C--     1           PPE,PSC,PTA,PWE,PSH,TSOIL(8)
cc      COMMON/FPMFG1/FGPM(15)
C--      COMMON/FPMFG1/itta,FGPM(15),ivers,ifrze
C--      COMMON/FRDSTFG/SMAX,PSISAT,BRT,SWLT,QUARTZ,STYPE,NUPL,
C--     +               NSAC,RTUZ,RTLZ,DZUP,DZLOW
cc      COMMON/FRZCNST/ FRST_FACT,CKSOIL,ZBOT
C--      COMMON/FRZCNST/ FRST_FACT,ZBOT      
      
      DO I=1,NSOIL
       STC(I)=FRZST(I)+T0
       ZSOIL(I)=FRZPAR(9+I)
      ENDDO

c 2/2012 error in estimation ITXT
c 2/2012 error      ITXT=FRZPAR(1)+0.01
      ITXT=FRZPAR(1)+0.5
      BRT = 15.9*CLAY(ITXT) + 2.91
      SMAX = -0.126*SAND(ITXT) + 0.489
      STYPE = 12

C  CALL SUBROUTINE TO RECALCULATE SAC-SMA STATES INTO FROZEN
C  GROUND MODEL STATES:
C  UPPER ZONE STORAGES
c      WRITE(*, *) 'FROST2_1: SACST(2) = ', SACST(2)
c      WRITE(*, *) 'FROST2_1: SACST_PRV(2) = ', SACST_PRV(2)
    
      DWT=FRZPAR(6)*(SACST(1)-SACST_PRV(1))
      DWF=FRZPAR(6)*(SACST(2)-SACST_PRV(2))
      NUP=2
cc      IEXIT=0
CVK 12/2005      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NUPL,ZSOIL,SMAX)
      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NUPL,ZSOIL,SMAX,FRZPAR(9))

C  LOWER ZONE STORAGES
      DWT=FRZPAR(7)*(SACST(3)-SACST_PRV(3))
      DWF=FRZPAR(7)*(SACST(4)+SACST(5)-SACST_PRV(4)-SACST_PRV(5))
      NUP=NUPL+1
CVK 12/2005      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NSAC,ZSOIL,SMAX)
      CALL SAC2FRZ1(DWT,DWF,SMC,SH2O,NUP,NSAC,ZSOIL,SMAX,FRZPAR(9))

c      WRITE(*,*) 'FROST2_1: ', (SMC(i), SH2O(i), i=1,NSAC)
C  CALCULATE SNOW DENSITY
C   SH - SNOW DEPTH IN CM
C   SR - SNOW DENSITY IN G/CM3
C   WE - SNOW WATER EQUIVALENT IN MM
      IF(WE .EQ. 0.) THEN
       SHX=0.
       SR=0.2
      ELSE
       SR=0.1*WE/SH
       SHX=0.01*SH
      ENDIF
C  CONVERT AIR TEMPERATURE INTO KALVIN UNITS
      TAIR=TA+T0    

C  FROZEN GROUND SIMULATION LOOP 
      DO IT=1,IDT
C  HRT1 ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN
       CALL HRT1(RHSTS,STC,SMC,SMAX,NSOIL,ZSOIL,TAIR,FRZPAR(2),
     +           FRZPAR(5),FRZPAR(8),SH2O,DTFRZ,BRT,SHX,SR,AESC,
     +           STYPE,SAND(ITXT),FRZPAR(3),IVERS,FRZPAR(4))
     
C     HSTEP ROUTINE CALCS/UPDATES SOIL TEMPS BASED ON RHSTS
       CALL HSTEP ( STCOUT,STC,RHSTS,DTFRZ,NSOIL)

C   DOUBLE CALLING HRT1 & HSTEP TO REDUCE NOISE. 
C   IT RUNS ONLY IF THERE IS AN ICE IN ONE SOIL LAYER AT LEAST
       IFRZ = 0
       DO I = 1,NSOIL
        XXS=SMC(I) - SH2O(I)
        IF( XXS .GT. 0.0 ) IFRZ = 1
       ENDDO
       IF ( IFRZ .EQ. 1 ) THEN     

C  SECOND CALL OF HRT1 AND HSTEP
        DO I = 1,NSOIL
         STC(I) = 0.5 * ( STCOUT(I) + STC(I) )
        ENDDO
        CALL HRT1(RHSTS,STC,SMC,SMAX,NSOIL,ZSOIL,TAIR,FRZPAR(2),
     +            FRZPAR(5),FRZPAR(8),SH2O,DTFRZ,BRT,SHX,SR,AESC,
     +            STYPE,SAND(ITXT),FRZPAR(3),IVERS,FRZPAR(4))
        
C     HSTEP ROUTINE CALCS/UPDATES SOIL TEMPS BASED ON RHSTS
        CALL HSTEP ( STC,STC,RHSTS,DTFRZ,NSOIL)
       ELSE

C  SKIP SECOND CALL OF HRT1 AND HSTEP
        DO I = 1,NSOIL
          STC(I) = STCOUT(I)
        ENDDO
       ENDIF

      ENDDO

C  STORE SOIL TEMPERATURE STATES IN CELSIUS
      ITFRZ=0
      DO I=1,NSOIL
       FRZST(I)=STC(I)-T0
       IF(FRZST(I) .LT. 0.) ITFRZ=ITFRZ+1
      ENDDO
      FRZDUP=0.
      FRZDBT=0.
      FROST=0.
      
C  RECALCULATE UNFROZEN MOISTURE STATES INTO SAC-SMA STORAGES 
      ISFRZ=0
      DO J=1,5
       IF(SACST(J) .NE. FRZST(J+5)) ISFRZ=ISFRZ+1
      ENDDO
      IF(ISFRZ .NE. 0 .OR. ITFRZ .NE. 0) THEN   
       CALL FRZ2SAC1(PX,ZSOIL,SMC,SH2O,NUPL,NSAC,FRZPAR(9),FRZPAR(6),
     +      FRZPAR(7),FRZST,SACST,FROST, HRAPX,HRAPY)
       IF(NSAC .NE. NSOIL) FROST=FROST-1000.*(SMC(NSOIL)-
     +    SH2O(NSOIL))*(ZSOIL(NSOIL-1)-ZSOIL(NSOIL))/FRST_FACT
     
C  FROST DEPTH CALCULATION
C--       FIND=FRZIND1(SMC,SH2O,TSOIL,ZSOIL,NSOIL,NUP,NSAC,RTUZ,RTLZ,
C--     +             FROST,SWLT,FRZ)
       CALL FRZIND1(SMC,SH2O,FRZST,ZSOIL,NSOIL,FRZPAR(9),FRZDUP,FRZDBT)
      ENDIF
       
C  REMEMBER SAC-SMA PREVIOUS STATES
      DO J=1,5
       SACST_PRV(J) = SACST(J)
      ENDDO 
C--      UZTWC0=UZTWC
C--      UZFWC0=UZFWC
C--      LZTWC0=LZTWC
C--      LZFSC0=LZFSC
C--      LZFPC0=LZFPC

CVK_02  NEW OPTION TO INTERPOLATE MODEL SOIL LAYER TEMP. INTO DESIRED LAYERS
C--      NMOD=NSOIL+1
CCC      DSMOD(1)=0.
C--      DSMOD(1)=-0.5*ZSOIL(1)
C--      TSTMP(1)=TSOIL(1)
C--      SWTMP(1)=SMC(2)
C--      SWHTMP(1)=SH2O(2)
C--      DSMOD(NMOD)=-ZBOT
C--      TSTMP(NMOD)=TBOT-T0
C--      SWTMP(NMOD)=SMAX
C--      SWHTMP(NMOD)=SMAX
C--      DO I=2,NMOD-1
C--       DSMOD(I)=-0.5*(ZSOIL(I-1)+ZSOIL(I))
C--       TSTMP(I)=TSOIL(I)
C--       SWTMP(I)=SMC(I)
C--       SWHTMP(I)=SH2O(I)
C--      ENDDO 
C--      CALL SOIL_INT1(TSTMP,NMOD,DSMOD,DSINT,NINT,TSINT)
C--      CALL SOIL_INT1(SWTMP,NMOD,DSMOD,DSINTW,NINTW,SWINT)
C--      CALL SOIL_INT1(SWHTMP,NMOD,DSMOD,DSINTW,NINTW,SWHINT)
C--      DO I=1,NINTW
C--       IF(I .EQ. 1) THEN
C--        SWINT(I)=SWINT(I)*DSINTW(I)*1000.
C--        SWHINT(I)=SWHINT(I)*DSINTW(I)*1000
C--       ELSE	
C--        SWINT(I)=SWINT(I)*(DSINTW(I)-DSINTW(I-1))*1000.
C--        SWHINT(I)=SWHINT(I)*(DSINTW(I)-DSINTW(I-1))*1000
C--       ENDIF
C--      ENDDO 	

      RETURN
      END
	  
	  