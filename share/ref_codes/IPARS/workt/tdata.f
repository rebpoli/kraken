C  TDATA.F - PROCESS TRANSIENT DATA
C          - SET DEFAULT EXTERNAL UNITS

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE GETTDAT (NERR)
C  SUBROUTINE GETTIMD (NERR)
C  SUBROUTINE GETTSCL (NERR)
C  SUBROUTINE SETEXU  (NERR)

C  CODE HISTORY:

C  JOHN WHEELER     12/2/95     ALPHA CODE
C  JOHN WHEELER      5/7/96     IMPLEMENT EXTERNAL UNITS
C  JOHN WHEELER     4/28/96     READ LINEAR SOLVER TRANSIENT DATA
C  JOHN WHEELER    11/10/98     REVISE EXTERNAL WELL RATE UNITS
C  RICK DEAN         3/3/01     SET UNITS FOR CONCENTRATION FOR COMPOSITIONAL
C                               SETUP WELL UNITS FOR COMPOSITIONAL
C  GURPREET SINGH   11/11/2015  ADDED 2PHASE IMPLICIT MFMFE
C  GERGINA PENCHEVA  4/22/2016  ADDED STOP-ON-ERROR STATEMENTS AND ADDITIONAL
C                               WARNINGS ON INPUT DATA (IN THE OUTPUT FILE)

C  NOTES:

C  1)  Subroutine GETTDAT should be called by all models; thus it
C      should contain no model specific input but may contain model
C      specific calls to other routines.

C*********************************************************************
      SUBROUTINE GETTDAT (NERR)
C*********************************************************************

C  Executive routine for extracting transient data from a keyword super
C  array.

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************

      INCLUDE 'control.h'

C  GET TIME CONTROL DATA

      CALL GETTIMD (NERR)
      IF (NERR.GT.0) STOP 'Error in GETTIMD'

C  GET FRAMEWORK DATA

      CALL GETTSCL (NERR)
      IF (NERR.GT.0) STOP 'Error in GETTSCL'

C  GET LINEAR SOLVER DATA

C        CALL LSORI (NERR,2)
C      CALL MULGRDI (NERR,2)
C      CALL SOLI (NERR,2)
C      CALL SOLGI (NERR,2)
C      CALL PCGI (NERR,2)
      IF (NERR.GT.0) STOP 'Error in solver intialization'

C  GET WELL DATA (INCLUDING MODEL SPECIFIC DATA)

      CALL IWELLS (NERR,2)
      IF (NERR.GT.0) STOP 'Error in IWELLS'
      IF (LEVERR.GT.2) RETURN

C  GET MODEL SPECIFIC TRANSIENT DATA

C      MODACT=15
C      IF (MODELON(15)) THEN
C         CALL ETDATA(NERR)
C         GOTO 1
C      ENDIF

C      MODACT=14
C      IF (MODELON(14)) THEN
C         CALL TRTDATA(NERR)
C         GOTO 1
C      ENDIF

C      MODACT=2
C      IF (MODELON(2)) CALL ITDATA(NERR)
C      MODACT=3
C      IF (MODELON(3)) CALL XTDATA(NERR)
C      MODACT=16
C      IF (MODELON(16)) CALL XTDATA(NERR)
C      MODACT=5
C      IF (MODELON(5)) CALL HTDATA(NERR)
C      MODACT=19
C      IF (MODELON(19)) CALL HTDATA(NERR)
C       MODACT=18
C       IF (MODELON(18)) CALL HTDATA(NERR)
C      MODACT=7
C      IF (MODELON(7)) CALL MTDATA(NERR)
      MODACT=13
      IF (MODELON(13)) CALL TTDATA(NERR)
C      MODACT=17
C      IF (MODELON(17)) CALL TTDATA(NERR)

    1 MODACT=0

C  CHECK FOR UNDEFINED DATA

      CALL UNDEF(NERR)
      IF ((NERR.GT.0).AND.LEVELC)
     & WRITE (NFOUT,*) 'TOTAL TRANSIENT DATA ERRORS =',NERR

      END
C*********************************************************************
      SUBROUTINE GETTIMD (NERR)
C*********************************************************************

C  EXTRACT TIME STEP AND OUTPUT CONTROL DATA FROM KEY-WORD INPUT

C  NERR = ERROR NUMBER STEPED BY 1 ON ERROR (INPUT & OUTPUT, INTEGER)

C*********************************************************************

      INCLUDE 'control.h'
      INCLUDE 'unitsex.h'
      INCLUDE 'restc.h'

      LOGICAL BK(15)
      CHARACTER*50 TITL

      IF (LEVELC) THEN
         TITL='TRANSIENT INPUT DATA'
         WRITE (NFOUT,*)
         CALL PRTTIT (TITL)
         WRITE (NFOUT,1) TIM*CVMTIME,EXTTIME
    1    FORMAT(/' TIME =',T45,F12.3,1X,A20)
         WRITE (NFOUT,10) NSTEP
   10    FORMAT(' TIME STEPS COMPLETED =',T45,I12)
      ENDIF

C  INPUT TIME STEP CONTROL

      CALL GETVAL('DELTIM[day] ',DELTIM,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,2) DELTIM*CVMTIME,EXTTIME
    2 FORMAT(' TARGET TIME STEP SIZE (DELTIM) =',T45,F12.3,1X,A20)

      CALL GETVAL('DTIMMUL ',DTIMMUL,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,3) DTIMMUL
    3 FORMAT(' TIME STEP MULTIPLIER (DTIMMUL) =',T45,F12.3)

      CALL GETVAL('DTIMMAX[day] ',DTIMMAX,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,4) DTIMMAX*CVMTIME,
     & EXTTIME
    4 FORMAT(' MAXIMUM TIME STEP (DTIMMAX) =',T45,F12.3,1X,A20)

      CALL GETVAL('DTIMMIN[day] ',DTIMMIN,'R8',0,0,0,0,NDUM,NERR)
      IF (NDUM.GT.0) THEN
         DTIMTOL=.1D0*DTIMMIN
         IF (LEVELC) WRITE (NFOUT,9) DTIMMAX*CVMTIME,EXTTIME
    9    FORMAT(' MINIMUM TIME STEP (DTIMMIN) =',T45,F12.3,1X,A20)
      ENDIF

C  INPUT PRINT CONTROL

      CALL GETVAL('TIMOUT[day] ',ACTTIM(3),'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,5) ACTTIM(3)*CVMTIME,
     & EXTTIME
    5 FORMAT(' NEXT SCHEDULED OUTPUT (TIMOUT) =',T45,F12.3,1X,A20)

      CALL GETVAL('DTIMOUT[day] ',DTIMOUT,'R8',0,0,0,0,NDUM,NERR)
      IF (DTIMOUT.LT.DTIMMIN) THEN
         DTIMOUT=DTIMMIN
         NDUM=1
      ENDIF

      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,6) DTIMOUT*CVMTIME,
     & EXTTIME
    6 FORMAT(' TIME BETWEEN OUTPUTS (DTIMOUT) =',T45,F12.3,1X,A20)

      CALL GETVAL('OUTLEVEL ',L,'I4',0,0,0,0,NDUM,NERR)
      IF (NDUM.GT.0) THEN
         LEVELA=.FALSE.
         LEVELB=.FALSE.
         LEVELC=.FALSE.
         IF (MYPRC.EQ.0) THEN
            LEVELC=.TRUE.
            IF (L.GT.1) LEVELB=.TRUE.
            IF (L.GT.2) LEVELA=.TRUE.
            IF (L.EQ.1) WRITE (NFOUT,17)
            IF (L.EQ.2) WRITE (NFOUT,18)
            IF (L.EQ.3) WRITE (NFOUT,19)
         ENDIF
   17    FORMAT(' USER OUTPUT LEVEL (OUTLEVEL) =',T56,'MINIMUM')
   18    FORMAT(' USER OUTPUT LEVEL (OUTLEVEL) =',T56,'STANDARD')
   19    FORMAT(' USER OUTPUT LEVEL (OUTLEVEL) =',T56,'MAXIMUM')
      ENDIF

      CALL GETVAL('DEBUGS ',LEVELD,'L4',0,0,0,0,NDUM,NERR)
      IF (NDUM.GT.0) THEN
         IF (LEVELD) THEN
            WRITE (NFOUT,26)
   26       FORMAT(' SINGLE PROCESSOR DEBUG (DEBUGS) =',T56,'YES')
         ELSE
            WRITE (NFOUT,27)
   27       FORMAT(' SINGLE PROCESSOR DEBUG (DEBUGS) =',T56,'NO')
         ENDIF
      ENDIF

      CALL GETVAL('DEBUGM ',LEVELE,'L4',0,0,0,0,NDUM,NERR)
      IF (NDUM.GT.0) THEN
         IF (LEVELE) THEN
            WRITE (NFOUT,28)
   28       FORMAT(' MULTI-PROCESSOR DEBUG (DEBUGM) =',T56,'YES')
         ELSE
            WRITE (NFOUT,29)
   29       FORMAT(' MULTI-PROCESSOR DEBUG (DEBUGM) =',T56,'NO')
         ENDIF
      ENDIF

      DO 11 I=1,15
   11 BK(I)=.FALSE.
      CALL GETVAL('BUGKEY ',BK,'FG',15,0,0,0,NDUM,NERR)
      DO 12 I=1,15
      IF (BK(I)) THEN
         BUGKEY(I)=.NOT.BUGKEY(I)
         IF (BUGKEY(I)) CALL OPENBUG()
      ENDIF
   12 CONTINUE

C  INPUT RESTART CONTROL

      CALL GETVAL('TIMRES[day] ',ACTTIM(4),'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,7) ACTTIM(4)*CVMTIME,
     & EXTTIME
    7 FORMAT(' NEXT RESTART OUTPUT (TIMRES) =',T45,F12.3,1X,A20)

      CALL GETVAL('DTIMRES[day] ',DTIMRES,'R8',0,0,0,0,NDUM,NERR)
      IF (DTIMRES.LT.DTIMMIN) THEN
         DTIMRES=DTIMMIN
         NDUM=1
      ENDIF
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,8) DTIMRES*CVMTIME,
     & EXTTIME
    8 FORMAT(' TIME BETWEEN RESTART OUTPUTS (DTIMRES) =',T45,F12.3,
     & 1X,A20)

      CALL GETVAL('FORMAT ',FORMOUT,'L4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,30) FORMOUT
   30 FORMAT(' RESTART OUTPUT FORMATTED (FORMAT) =',T46,L8)

C  INPUT VISUALIZATION OUTPUT CONTROL

      CALL GETVAL('VISALL ',VISALL,'FG',0,0,0,0,NDUM,NERR)
      CALL GETVAL('VISNEWT ',VISNEWT,'FG',0,0,0,0,NDUM,NERR)

      IF (VISNEWT) THEN
         IF (LEVELC) WRITE (NFOUT,*)
     &        'VISALL: VIS.output at every NEWTON STEP. '

      ELSEIF (VISALL) THEN
         IF (LEVELC) WRITE (NFOUT,*)
     &        'VISALL: VIS.output at every TIME STEP. '

      ELSE
         CALL GETVAL('VISOUT[day] ',ACTTIM(6),'R8',0,0,0,0,NDUM,NERR)
         IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,59) ACTTIM(6)*CVMTIME,
     &        EXTTIME
 59   FORMAT(' NEXT SCHEDULED VISUALIZATION OUTPUT =',T45,F12.3,1X,A20)

         CALL GETVAL('DVISOUT[day] ',DVISOUT,'R8',0,0,0,0,NDUM,NERR)
         IF (DVISOUT.LT.DTIMMIN) THEN
            DVISOUT=DTIMMIN
            NDUM=1
         ENDIF
      ENDIF

      call GET_VISPARAMS(nerr)

      IF (LEVELC.AND.(NDUM.GT.0)) WRITE (NFOUT,510) DVISOUT*CVMTIME,
     & EXTTIME
c23456
  510  FORMAT(' TIME BETWEEN SCHEDULED VISUALIZATION OUTPUTS =',
     &  T45,F12.3,1X,A20)

      END
C*********************************************************************
      SUBROUTINE GETTSCL (NERR)
C*********************************************************************

C  EXTRACT MISCELLENOUS TRANSIENT DATA FROM KEY-WORD INPUT

C  NERR = ERROR NUMBER STEPED BY 1 ON ERROR (INPUT & OUTPUT, INTEGER)

C*********************************************************************

C  GET REVISED ARRAY OUTPUT INDEXES

      CALL GETGAPI (.FALSE.,NERR)

      END
C*********************************************************************
      SUBROUTINE SETEXU (NERR)
C*********************************************************************

C  INITIALIZES CONVERSION BETWEEN INTERNAL AND EXTERNAL UNITS

C  NERR = ERROR NUMBER STEPED BY 1 ON ERROR (INPUT & OUTPUT, INTEGER)

C Value External = (Value Internal) * CVMxxxx + CVAxxxx

C External External     Internal Internal    Definition
C   Name   Initial        Name   Initial
C          Default               Default

C EXTMASS  [lb]         INTMASS  [lb]        Mass
C EXTDIST  [ft]         INTDIST  [ft]        Distance
C EXTVOLL  [bbl]        INTVOLL  [cu-ft]     Volume
C EXTTIME  [day]        INTTIME  [day]       Time
C EXTTEMP  [F]          INTTEMP  [F]         Temperature
C EXTPRES  [psi]        INTPRES  [psi]       Pressure
C EXTPERM  [md]         INTPERM  [md]        Permeability
C EXTWELL  [stb/day]    INTWELL  [lb/day]    Liquid mass rate
C EXTWELG  [mscf/day]   INTWELG  [lb/day]    Gas mass rate
C EXTWELX  [lbM/day]    INTWELX  [lbM/day]   Component mass rate
C EXTWELXC [lbM]                             Component mass
C EXTVISC  [cp]         INTVISC  [cp]        Viscosity
C EXTDENS  [lb/cu-ft]   INTDENS  [lb/cu-ft]  Density
C EXTCONC  [M/cu-ft]    INTCONC  [M/cu-ft]   Molar concentration
C EXTCOMP  [/psi]       INTCOMP  [/psi]      Compressability

C*********************************************************************
      INCLUDE 'control.h'
      INCLUDE 'unitsex.h'

      IF (MYPRC.EQ.0) THEN
         WRITE (NFOUT,*)
         WRITE (NFOUT,*) 'DEFAULT EXTERNAL UNITS:'
      ENDIF

C  SET DEFAULT INTERNAL AND EXTERNAL UNITS

      EXTMASS='[lb]'
      INTMASS='[lb]'
      EXTDIST='[ft]'
      INTDIST='[ft]'
      EXTVOLL='[bbl]'
      INTVOLL='[cu-ft]'
      EXTTIME='[day]'
      INTTIME='[day]'
      EXTTEMP='[F]'
      INTTEMP='[F]'
      EXTPRES='[psi]'
      INTPRES='[psi]'
      EXTPERM='[md]'
      INTPERM='[md]'
      EXTWELL='[stb/day]'
      EXTWELLC='[stb]'
      INTWELL='[lb/day]'
      EXTWELG='[mscf/day]'
      EXTWELGC='[mscf]'
      INTWELG='[lb/day]'
      EXTVISC='[cp]'
      INTVISC='[cp]'
      EXTDENS='[lb/cu-ft]'
      INTDENS='[lb/cu-ft]'
      EXTCONC='[M/cu-ft]'
      INTCONC='[M/cu-ft]'
      EXTCOMP='[/psi]'
      INTCOMP='[/psi]'

C      EXTMASS='[lbM]'
C      INTMASS='[lbM]'
C      EXTCONC='[lbM/cu-ft]'
C      INTCONC='[lbM/cu-ft]'
C      INTWELL='[stb/day]'
C      INTWELG='[mscf/day]'
C      INTWELX='[lbM/day]'
C      EXTWELX='[lbM/day]'
C      EXTWELXC='[lbM]'


C      EXTMASS='[lbM]'
C      INTMASS='[lbM]'
C      EXTCONC='[lbM/cu-ft]'
C      INTCONC='[lbM/cu-ft]'
C      INTWELL='[stb/day]'
C      INTWELG='[mscf/day]'
C      INTWELX='[lbM/day]'
C      EXTWELX='[lbM/day]'
C      EXTWELXC='[lbM]'

cgp
C      EXTCONC='[lbM/cu-ft]'
C      INTCONC='[lbM/cu-ft]'

C  GET STOCK TANK DENSITIES AND SET CONVERSION FACTORS

      STDENW=62.34D0
      CALL GETVAL('STDENW[lb/cu-ft] ',STDENW,'R8',0,0,0,0,NDUM,NERR)
C     2.546732D0 kg cu-ft / lb bbl
      CALL SETCVF('stb          ',STDENW*2.546732D0)
      CALL SETCVF('stbw         ',STDENW*2.546732D0)

      STDENO=56.D0
      CALL GETVAL('STDENO[lb/cu-ft] ',STDENO,'R8',0,0,0,0,NDUM,NERR)
      CALL SETCVF('stbo         ',STDENO*2.546732D0)

      STDENG=.04228D0
      CALL GETVAL('STDENG[lb/cu-ft] ',STDENG,'R8',0,0,0,0,NDUM,NERR)
C     .4535924 kg / lb
      CALL SETCVF('scf          ',STDENG*.4535924D0)
      CALL SETCVF('mscf         ',STDENG*453.5924D0)
C     16.018466 kg cu-ft/ lb cu-m
      CALL SETCVF('scm          ',STDENG*16.018466D0)

C  GET EXTERNAL UNITS

      CALL EXTUIN('EXTMASS ',EXTMASS,INTMASS,CVMMASS,CVAMASS,NERR,
     & 'MASS                ')

      CALL EXTUIN('EXTDIST ',EXTDIST,INTDIST,CVMDIST,CVADIST,NERR,
     & 'DISTANCE            ')

      CALL EXTUIN('EXTVOLL ',EXTVOLL,INTVOLL,CVMVOLL,CVAVOLL,NERR,
     & 'VOLUME              ')

      CALL EXTUIN('EXTTIME ',EXTTIME,INTTIME,CVMTIME,CVATIME,NERR,
     & 'TIME                ')

      CALL EXTUIN('EXTTEMP ',EXTTEMP,INTTEMP,CVMTEMP,CVATEMP,NERR,
     & 'TEMPERATURE         ')

      CALL EXTUIN('EXTPRES ',EXTPRES,INTPRES,CVMPRES,CVAPRES,NERR,
     & 'PRESSURE            ')

      CALL EXTUIN('EXTPERM ',EXTPERM,INTPERM,CVMPERM,CVAPERM,NERR,
     & 'PERMEABILITY        ')

      CALL EXTUIN('EXTWELL ',EXTWELL,INTWELL,CVMWELL,CVAWELL,NERR,
     & 'LIQUID MASS RATE    ')
      CALL INSTR(3,'stb',20,EXTWELL,N)
      STBEXT=.FALSE.
      IF (N.GT.0) STBEXT=.TRUE.

      CALL EXTUIN('EXTWELG ',EXTWELG,INTWELG,CVMWELG,CVAWELG,NERR,
     & 'GAS MASS RATE       ')
      CALL INSTR(3,'scf',20,EXTWELG,N)
      SCFEXT=.FALSE.
      IF (N.GT.0) SCFEXT=.TRUE.
      CALL INSTR(3,'scm',20,EXTWELG,N)
      IF (N.GT.0) SCFEXT=.TRUE.

      CALL EXTUIN('EXTVISC ',EXTVISC,INTVISC,CVMVISC,CVAVISC,NERR,
     & 'VISCOSITY           ')

      CALL EXTUIN('EXTDENS ',EXTDENS,INTDENS,CVMDENS,CVADENS,NERR,
     & 'DENSITY             ')

      CALL EXTUIN('EXTCONC ',EXTCONC,INTCONC,CVMCONC,CVACONC,NERR,
     & 'MOLAR CONCENTRATION ')

      CALL EXTUIN('EXTCOMP ',EXTCOMP,INTCOMP,CVMCOMP,CVACOMP,NERR,
     & 'COMPRESSIBILITY     ')

      IF (LEVELC) THEN

         WRITE (NFOUT,1) STDENO*CVMDENS,EXTDENS
    1    FORMAT(/' STOCK TANK OIL DENSITY (STDENO)',T48,F12.3,1X,A20)
         WRITE (NFOUT,2) STDENW*CVMDENS,EXTDENS
    2    FORMAT(' STOCK TANK WATER DENSITY (STDENW)',T48,F12.3,1X,A20)
         WRITE (NFOUT,3) STDENG*CVMDENS,EXTDENS
    3    FORMAT(' STOCK TANK GAS DENSITY (STDENG)',T48,F12.5,1X,A20)

      ENDIF

      END
C*********************************************************************
      SUBROUTINE EXTUIN (KNAM,EXT,INTU,CVM,CVA,NERR,DISC)
C*********************************************************************

C  SERVICE ROUTINE FOR SETEXU()

C  KNAM = EXTERNAL UNITS KEYWORD (INPUT, CHARACTER*8)

C  EXT  = EXTERNAL UNITS STRING (INPUT & OUTPUT, CHARACTER*20)

C  INTU = INTERNAL UNITS STRING (INPUT, CHARACTER*20)

C  CVM  = C.F. MULTIPLIER FROM INTERNAL TO EXTERNAL UNITS (OUTPUT, REAL*8)

C  CVA  = C.F. ADDITION FROM INTERNAL TO EXTERNAL UNITS (OUTPUT, REAL*8)

C  NERR = ERROR NUMBER STEPED BY 1 ON ERROR (INPUT & OUTPUT, INTEGER)

C  DISC = DISCRIPTION STRING (INPUT, CHARACTER*20)
C*********************************************************************
      INCLUDE 'control.h'

      REAL*8 CVM,CVA
      CHARACTER*8 KNAM
      CHARACTER*20 EXT,INTU,DISC
      CHARACTER*50 ERRMSG

      CALL GETVALS(KNAM,EXT,'CS',0,0,0,20,NDUM,NERR)
      IF (MYPRC.EQ.0) WRITE (NFOUT,1) DISC,KNAM,EXT
    1 FORMAT (4X,A20,'(',A8,')',T56,A20)
      CALL CONVRT (INTU,EXT,CVM,CVA,KODRET,ERRMSG)
      IF (KODRET.GT.0) THEN
         NERR=NERR+1
         IF (MYPRC.EQ.0) WRITE (NFOUT,2) KODRET,ERRMSG,INTU,EXT
    2    FORMAT(/' ERROR #',I4,' - ',A50/' CONVERTING FROM ',A20,
     &   ' TO ',A20)
      ENDIF
      END
