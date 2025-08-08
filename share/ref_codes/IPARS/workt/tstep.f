C  TSTEP.F - MAKE ONE TIME STEP WITH THE IMPLICIT SINGLE PHASE FLOW MODEL

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE TSTEP1   (NERR)
C  SUBROUTINE TSTEP2   (KONVG,NERR)
C  SUBROUTINE TSTEP3   (KONVG,NERR)

C  CODE HISTORY:

C  BAHAREH MOMKEN 3/16/99  hydrology(IMPES) gstep.df is used as template
C  JOHN WHEELER   4/03/99  IMPLICIT SINGLE PHASE MODEL
C  JOHN WHEELER   6/22/99  SPLIT TSTEP TO 3 ROUTINES FOR MULTIMODEL

C ***********************************************************************
        SUBROUTINE TSTEP1 (NERR)
C************************************************************************

C  Implicit single-phase model exective routine to start a time step.
C  1.  Compute Jacobian and residual for Newtonian iteration.
C  2.  Start inter-block calculations.
C
C  NERR = ERROR KEY STEPPED BY ONE FOR EACH ERROR
C         (INPUT AND OUTPUT, INTEGER )
C*********************************************************************
      IMPLICIT NONE
C    INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'wells.h'
      INCLUDE 'tarydat.h'
      INCLUDE 'tbaldat.h'
      INCLUDE 'terrcalc.h'

      LOGICAL ONCEONLY
      DATA ONCEONLY /.TRUE./

      INTEGER NERR
      INTEGER IPROP(8), ITRAN(9), IWELL(6), IWELSUM(4)

      DATA IPROP   / 8*0 /
      DATA ITRAN   / 9*0 /
      DATA IWELL   / 6*0 /
      DATA IWELSUM / 4*0 /

c list all routines to be called by callwork:

      EXTERNAL TPROP, TTRAN, TWELL, TWELSUMS
C      EXTERNAL TSPRB3

c  define work routine arguments

      IF(ONCEONLY) THEN
         ONCEONLY = .FALSE.

         IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)'PROC',MYPRC,
     &    ' ENTERING SUBROUTINE TSTEP1, OLD TAG =',MSGTAG(13+1)

         IPROP(1) = 7
         IPROP(2) = N_FLDEN
         IPROP(3) = N_PRES
         IPROP(4) = N_FLDENN
         IPROP(5) = N_PRESN
         IPROP(6) = N_COF
         IPROP(7) = N_RESID
         IPROP(8) = N_POR

         ITRAN(1) = 8
         ITRAN(2) = N_TCOFX
         ITRAN(3) = N_TCOFY
         ITRAN(4) = N_TCOFZ
         ITRAN(5) = N_DEPTH
         ITRAN(6) = N_FLDEN
         ITRAN(7) = N_PRES
         ITRAN(8) = N_COF
         ITRAN(9) = N_RESID

         IWELL(1) = 5
         IWELL(2) = N_DEPTH
         IWELL(3) = N_PRES
         IWELL(4) = N_FLDEN
         IWELL(5) = N_COF
         IWELL(6) = N_RESID

         IWELSUM(1) = 3
         IWELSUM(2) = N_DEPTH
         IWELSUM(3) = N_PRES
         IWELSUM(4) = N_FLDEN

      ENDIF

C  EVALUATE PROPERTIES ON DIRICHLET BOUNDARIES

      CALL TIMON(18)
      CALL TBDPROP()
      CALL TIMOFF(18)

C  START NEWTONIAN ITERATION LOOP
C  EVALUATE DENSITY AND ACCUMULATION TERMS

      CALL TIMON(21)
      CURRENT=0.D0
      CALL CALLWORK(TPROP,IPROP)
      CALL TIMOFF(21)
      IF((NHISUSE == 0).AND.(NSTEP < 1)) GO TO 1

      IF (NEWT.EQ.1.AND.TIM.EQ.0.D0) BALANCE(1,MODACT,4)=CURRENT

C  EXCHANGE PHYSICAL PARAMETERS WITH NEIGHBORING PROCESSORS

      CALL TIMON(22)
      CALL UPDATE(N_PRES,1)
      CALL UPDATE(N_FLDEN,1)
      CALL TIMOFF(22)

C  EVALUATE TRANSPORT IN COEFFICIENTS AND RESIDUALS

      CALL TIMON(21)
      CALL CALLWORK(TTRAN,ITRAN)
C         CALL CALLTSPRB3(TSPRB3,ITRAN)
      CALL TIMOFF(21)

C  EVALUATE TRANSPORT ACROSS DIRICHLET BOUNDARIES

      CALL TIMON(18)
      CALL TBDTRAN ()
      CALL TIMOFF(18)

C bag8 - source function for manufactured solution
      IF (ITEST.GT.0) CALL TSOURCE()

C evaluate well conditions and add well conditions to the
C pressure matrix and residual

    1 CONTINUE
      CALL TIMON(10)
      CALL CLEARWELLS()
      CALL CALLWORK(TWELSUMS,IWELSUM)
      CALL PARALLWELLS()
      CALL CALLWORK(TWELL,IWELL)
      CALL TIMOFF(10)

      IF((NHISUSE == 0).AND.(NSTEP < 1)) RETURN

C  SEND BLOCK INTERFACE DATA

      CALL TDUALS()

      END

C ***********************************************************************
        SUBROUTINE TSTEP2 (KONVG,NERR)
C************************************************************************
C  Implicit single-phase model executive routine to continue a time step
C  1.  Complete inter-block calculations
C  2.  Check Newtonian convergence
C
C  NERR = ERROR KEY STEPPED BY ONE FOR EACH ERROR
C         (INPUT AND OUTPUT, INTEGER )

C  KONVG = CONVERGENCE FLAG (OUTPUT, INTEGER)
C        = 1 ==> CONVERGED
C        = 2 ==> CONTINUE ITERATION
C        = 3 ==> FAILED
C        = 4 ==> CUT TIME STEP (SET BY FRAMEWORK ONLY)
C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'wells.h'
      INCLUDE 'layout.h'
      INCLUDE 'tarydat.h'
      INCLUDE 'tbaldat.h'

      INTEGER KONVG,NERR,IW,I
      EXTERNAL TBUGOUT, TMAXRESID

      LOGICAL ONCEONLY
      DATA ONCEONLY /.TRUE./
      INTEGER IBUGOUT(9), IMXRESID(2), NUMBUG6
      DATA NUMBUG6 /0/, IBUGOUT / 9*0 /, IMXRESID/ 2*0 /

      IF((NHISUSE == 0).AND.(NSTEP < 1)) RETURN

c  define work routine arguments

      IF(ONCEONLY) THEN
         ONCEONLY = .FALSE.

         IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)'PROC',MYPRC,
     &    ' ENTERING SUBROUTINE TSTEP2, OLD TAG =',MSGTAG(13+1)

         IBUGOUT(1) = 8
         IBUGOUT(2) = N_POR
         IBUGOUT(3) = N_PRES
         IBUGOUT(4) = N_FLDEN
         IBUGOUT(5) = N_COF
         IBUGOUT(6) = N_RESID
         IBUGOUT(7) = N_TCOFX
         IBUGOUT(8) = N_TCOFY
         IBUGOUT(9) = N_TCOFZ

         IMXRESID(1) = 1
         IMXRESID(2) = N_RESID

      ENDIF

C  RECEIVE BLOCK INTERFACE DATA

      FLUXOM=0.D0
      CALL TDUALR()

C  DEBUG OUTPUT

      IF (BUGKEY(6).AND.(NEWT.LT.3).AND.(NUMBUG6.LT.7)) THEN
         CALL CALLWORK(TBUGOUT,IBUGOUT)
         NUMBUG6=NUMBUG6+1
         IF (LEVELC) THEN
            DO 6 IW=1,NUMWEL
            IF (MODWEL(IW).EQ.MODACT
C     &     .OR. MODACT.EQ.FLOWMODEL
     &         )
     &         WRITE (NFBUG,7) IW,(WELIPC(I,IW),I=1,3)
    6       CONTINUE
         ENDIF
    7    FORMAT(' WELL',I3,' QWI',G11.4,' QWP',G11.4,' BHP',F8.2)
      ENDIF

C  CHECK NEWTONIAN CONVERGENCE

      RESIDT=0.D0
      RESIDM=0.D0
      RESID2=0.D0
      KONVG=2
      IF (NEWT.EQ.1) THEN
! bag8
         IF (DEALII) THEN
           CVTOL2=CVTOL
           CVTOL1=CVTOL
         ELSE
           CVTOL2=BALANCE(1,MODACT,4)*CVTOL
           CVTOL1=5.D0*CVTOL2/NEMOD(13)
         ENDIF
      ELSE
         CALL CALLWORK(TMAXRESID,IMXRESID)

         CALL MAXIT(1,RESIDM)
         CALL SPREAD8(1,RESIDM)
         CALL SUMIT(1,RESIDT)
         CALL SPREAD8(1,RESIDT)
         CALL SUMIT(1,RESID2)
         CALL SPREAD8(1,RESID2)
         RESID2=SQRT(RESID2)

         IF (NEWT.LE.MAXITS) THEN
            IF (RESIDM.LT.CVTOL1.AND.ABS(RESIDT).LT.CVTOL2) KONVG=1
         ELSE
            IF (RESIDM.LT.5.D0*CVTOL1.AND.
     &         ABS(RESIDT).LT.5.D0*CVTOL2) THEN
               KONVG=1
            ELSE
               KONVG=3
            ENDIF
         ENDIF
      ENDIF

      IF (BUGKEY(7)) THEN
         IF (LEVELC) WRITE (NFBUG,14) NEWT-1,RESIDM,RESIDT
         IF (MYPRC.EQ.0) WRITE (*,14) NEWT-1,RESIDM,RESIDT
   14    FORMAT(' NEWT =',I3,' Rmax =',G10.3,' Rtot =',G10.3)
      ENDIF

! bag8
      IF (MYPRC.EQ.0) THEN
        WRITE(*,10,advance="no") NEWT-1,NSTEP
        IF (RESIDM.LT.CVTOL1) THEN
          WRITE(*,30,advance="no") RESIDM,CVTOL1
        ELSE
          WRITE(*,40,advance="no") RESIDM,CVTOL1
        ENDIF
        IF (RESIDT.LT.CVTOL2) THEN
          WRITE(*,50) RESIDT,CVTOL2
        ELSE
          WRITE(*,60) RESIDT,CVTOL2
        ENDIF
      ENDIF
  10  FORMAT(' NEWT ',I4,' NSTEP ',I6,' NL-RESID')
  30  FORMAT(1P,' (MAX (',E9.2,' < ',E9.2,')')
  40  FORMAT(1P,' (MAX (',E9.2,' > ',E9.2,')')
  50  FORMAT(1P,' (TOT (',E9.2,' < ',E9.2,')')
  60  FORMAT(1P,' (TOT (',E9.2,' > ',E9.2,')')

      END

C ***********************************************************************
        SUBROUTINE TSTEP3 (KONVG,NERR)
C************************************************************************

C  Implicit single-phase model exective routine to complete a time step
C    1.  Update unknowns for next Newtonian iteration (KONVG = 2)
C    2.  Wrap up time step (KONVG = 1)
C    3.  Restart time step (KONVG = 4)
C    4.  Abort time step  (KONVG = 3)

C  KONVG = CONVERGENCE FLAG (INPUT, INTEGER)
C        = 1 ==> CONVERGED
C        = 2 ==> CONTINUE ITERATION
C        = 3 ==> FAILED
C        = 4 ==> CUT TIME STEP

C  NERR = ERROR KEY STEPPED BY ONE FOR EACH ERROR
C         (INPUT AND OUTPUT, INTEGER )
C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'wells.h'
      INCLUDE 'tarydat.h'
      INCLUDE 'tbaldat.h'
      INCLUDE 'tfluidsc.h'
      INCLUDE 'terrcalc.h'

      INTEGER KONVG,NERR
      EXTERNAL TUPPRES, CPYARYR8

      LOGICAL ONCEONLY
      DATA ONCEONLY /.TRUE./
      INTEGER ICOPY(3), IUPPRES(3)
      DATA ICOPY/ 3*0 /, IUPPRES / 3*0 /

c  define work routine arguments

      IF(ONCEONLY) THEN
         ONCEONLY = .FALSE.

         IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)'PROC',MYPRC,
     &    ' ENTERING SUBROUTINE TSTEP3, OLD TAG =',MSGTAG(13+1)

         ICOPY(1)=2
         ICOPY(2)=N_PRESN
         ICOPY(3)=N_PRES

         IUPPRES(1) = 2
         IUPPRES(2) = N_DUNK
         IUPPRES(3) = N_PRES

      ENDIF

      IF((NHISUSE == 0).AND.(NSTEP < 1)) GO TO 1

C  BRANCH ON CONVERGENCE STATE

      GO TO (1,2,1,4),KONVG

C  CONTINUE NEWTONIAN ITERATION

    2 IF (BUGKEY(6).AND.NEWT.LT.3) THEN
         TITU='CHANGE IN PRESSURE FOR FAULT BLOCK'
         CALL GEAOUT(N_DUNK,1,1)
      ENDIF

      CALL TIMON(20)
      CALL CALLWORK(TUPPRES, IUPPRES)
      CALL TIMOFF(20)

      RETURN

C  NEWTONIAN ITERATION CONVERGED OR FAILED
C  UPDATE WELL DATA AND BALANCES AT END OF TIME STEP

    1 CALL TWELLOUTPUT()
      IF((NHISUSE == 0).AND.(NSTEP < 1)) RETURN

      CALL TIMON(18)
      CALL TBDBAL ()
      CALL TIMOFF(18)

      BALANCE(1,MODACT,1)=CURRENT
      BALANCE(1,MODACT,2)=WELIS
      BALANCE(1,MODACT,3)=FLUXOM
      BALANCE(1,MODACT,7)=FLITNP

C  COMPUTE (POST-PROCESS) FACE-CENTERED VELOCITIES

      CALL TVEL_EX()

C  OUTPUT DIRICHLET BOUNDARY DATA

      CALL BND_OUT()

C  bag8 - ERROR CALCULATIONS

      IF (ITEST.GT.0) CALL TERRCALC()

      RETURN

C  TIME STEP CUT

    4 CALL CALLWORK(CPYARYR8,ICOPY)

      END
