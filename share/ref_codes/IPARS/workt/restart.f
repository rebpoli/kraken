C  RESTART.F - RESTART UTILITIES

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE RESOUT  (NERR)
C  SUBROUTINE RESIN   (NERR)
C  SUBROUTINE RGEAOUT (NARR,NEXTRA,NERR)
C  SUBROUTINE RSR4    (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,ARY)
C  SUBROUTINE RSR8    (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,ARY)
C  SUBROUTINE RSI4    (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,IARY)
C  SUBROUTINE RSL4    (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,LARY)
C  SUBROUTINE RGEAIN  (NERR)
C  SUBROUTINE RSINR4  (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,ARY)
C  SUBROUTINE RSINR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,ARY)
C  SUBROUTINE RSINI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,IARY)
C  SUBROUTINE RSINL4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK,IARY)

C  CODE HISTORY:

C  JOHN WHEELER      2/30/98    ALPHA CODE
C  JOHN WHEELER      4/20/99    ADD GENERIC MULTI-MODEL CAPABILITY
C  SUNIL G THOMAS    -/--/09    CHANGED MODACT TO CURRENT_MODEL AT
C                               TAG DEFINITION IN RGEAOUT(IN) TO KEEP
C                               PARALLEL PERFORMANCE IN FLOW COUPLED
C                               TO TRCHEM AND TO ACCOUNT FOR CHANGE
C                               IN PUTIL.DF (SEE COMMENT THERE). THE
C                               GENERIC MULTI-MODEL CAPABILITY CAN
C                               STILL BE PRESERVED AS ALLUDED THEREIN.

C  NOTES:

C  1.  Restart records will optionally be written in ASCII so that they are
C      portable.

C  2.  Restart records will be output between time steps before reading
C      transient data due at the restart time (if any).

C  3.  Restart will be accompished as follows:

C      A.  Initial data will be read and processed as normally done at
C          at time zero.

C      B.  Transient data up to but excluding the time of the restart will
C          be read and processed normally.

C      C.  Restart data will be read and processed.

C      D.  Transient data at the time of the restart (if any) will be read
C          and processed normally.

C  4.  Each record set will consist of a header record followed by data
C      records (optional)

C  5.  The format of header records will be (6I9) where the 1st integer
C      is the record set ID and the remaining 5 values are utility integers
C      which may be used to store restart data.  The 2nd integer in a
C      model header record must be the model number.

C  6.  Record set IDs are defined as follows

C      A.  1 ==> Framework control data

C      B.  2 ==> Last restart record

C      C.  3 ==> Start grid-element array data

C      E.  4 ==> Well data

C      F.  5 ==> Physical model data (other than grid-element array data)

C      G.  6 ==> Linear solver data

C  7.  The following field types are recommended but not required

C      A.  INTEGER      I9

C      B.  REAL         G15.8

C      C.  DOUBLE       G23.16

C      D.  STRING       I5              (STRING LENGTH)
C                       CHARACTER*1

C*********************************************************************
      SUBROUTINE RESOUT (NERR)
C*********************************************************************

C  ROUTINE DIRECTS OUTPUT OF A RESTART FILE

C  NERR = ERROR COUNT STEPPED BY ONE ON ERROR (INPUT AND OUTPUT, INTEGER)

C*********************************************************************
      INCLUDE 'control.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'

      CHARACTER*40 FN
      CHARACTER*10 EX
      CHARACTER*1 FN1(40),DOT,BLK,EX1(10),SPC,SUN,C
      EQUIVALENCE (FN,FN1),(EX,EX1)
      DATA DOT/'.'/,BLK/' '/

      SPC=CHAR(92)
      SUN=CHAR(47)

C  FORM RESTART FILENAME AND OPEN FILE (PROCESSOR 0)

      IF (MYPRC.NE.0) GO TO 15

C     PUT TIME STEP NUMBER IN EX

      WRITE (EX,1) NSTEP
    1 FORMAT(I10)
      M1=0
      DO 2 I=1,10
      IF (EX1(I).EQ.BLK) THEN
         IF (M1.GT.0) GO TO 3
      ELSE
         IF (M1.EQ.0) M1=I
      ENDIF
    2 M2=I
    3 M=M2-M1+1

C     COPY FILE NAME AND PICK OFF DIRECTORY NAME (IF ANY)

      L=0
      LS=0
      FN=' '
      DO 4 I=1,40
      C=RSFOUT1(I)
      IF (C.EQ.BLK) THEN
         IF (L.EQ.0) GO TO 4
         GO TO 5
      ENDIF
      L=L+1
      FN1(L)=C
      IF (C.EQ.SPC.OR.C.EQ.SUN) LS=L
    4 CONTINUE

    5 IF (L.EQ.0) THEN
         FN='R'
         L=1
      ENDIF

C     ADD EXTENSION AND OPEN FILE

      ML=LS+12-M
      IF (L.GT.ML) L=ML
      ML=L+1
      IF (M.GT.3) ML=ML-M+3
      IF (ML.LE.LS) THEN
         FN='R'
         L=1
         ML=2
      ENDIF
      FN1(ML)=DOT
      DO 6 I=1,M
    6 FN1(ML+I)=EX1(M1+I-1)

      IF (FORMOUT) THEN
         OPEN (NFROUT,FILE=FN,STATUS='UNKNOWN',ERR=13)
      ELSE
         OPEN (NFROUT,FILE=FN,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &    ERR=13)
      ENDIF

C  OUTPUT RESTART TIME AT TOP OF FILE

      IF (FORMOUT) THEN
         WRITE (NFROUT,9) TIM
      ELSE
         WRITE (NFROUT) TIM
      ENDIF

C  OUTPUT CONTROL RESTART RECORD SET (L=1)

      L=1
      WRITE (NFOUT,7) FN,NSTEP,TIM
      WRITE (*,7) FN,NSTEP,TIM
    7 FORMAT(' WRITING RESTART FILE ',A40/
     & '      AT STEP',I6,' AND TIME',F12.3)
      I2=3
      J2=19
      IF (FORMOUT) THEN
         WRITE (NFROUT,8) L,NSTEP,I2,J2,L,L
    8    FORMAT(6I9)
         WRITE (NFROUT,9) DELTIM,(ACTTIM(I),I=3,NACTTIM)
    9    FORMAT(5G23.16)
         WRITE (NFROUT,9) (((BALANCE(I,J,K),I=1,I2),J=1,J2),K=1,8)
      ELSE
         WRITE (NFROUT) L,NSTEP,I2,J2,L,L
         WRITE (NFROUT) DELTIM,(ACTTIM(I),I=3,NACTTIM)
         WRITE (NFROUT) (((BALANCE(I,J,K),I=1,I2),J=1,J2),K=1,8)
      ENDIF

C  OUTPUT FRAMEWORK RESTART ARRAYS (IF ANY) (L=3)

   15 CONTINUE

C  OUTPUT WELL RESTART DATA (IF ANY) (L=4)

C  OUTPUT LINEAR SOLVER DATA (IF ANY) (L=6)

C      CALL LSORRO(NERR)

C  OUTPUT PHYSICAL MODEL RESTART DATA (L=5)

C       MODACT=1
C       IF (MODELON(1)) CALL PRESTO(NERR)
C       MODACT=2
C       IF (MODELON(2)) CALL IRESTO(NERR)
C         MODACT=3
C         IF (MODELON(3)) CALL XRESTO(NERR)
C       MODACT=16
C       IF (MODELON(16)) CALL XRESTO(NERR)
C         MODACT=4
C         IF (MODELON(4)) CALL CRESTO(NERR)
C       MODACT=5
C       IF (MODELON(5)) CALL HRESTO(NERR)
C       MODACT=19
C       IF (MODELON(19)) CALL HRESTO(NERR)
C       MODACT=18
C       IF (MODELON(18)) CALL HRESTO(NERR)
C       MODACT=6
C       IF (MODELON(6)) CALL GRESTO(NERR)
C       MODACT=7
C       IF (MODELON(7)) CALL MRESTO(NERR)
C      MODACT=9
C      IF (MODELON(9)) CALL SRESTO(NERR)
      MODACT=13
      IF (MODELON(13)) CALL TRESTO(NERR)
C      MODACT=17
C      IF (MODELON(17)) CALL TRESTO(NERR)

              MODACT=0

C  OUTPUT TERMINAL RESTART RECORD AND EXIT (L=2)

      L=2
      IF (MYPRC.EQ.0) THEN
         IF (FORMOUT) THEN
            WRITE (NFROUT,8) L,L,L,L,L,L
         ELSE
            WRITE (NFROUT) L,L,L,L,L,L
         ENDIF
         CLOSE (NFROUT)
      ENDIF

      CALL WAITALL()
      RETURN

   13 WRITE (NFOUT,14)
   14 FORMAT (/' ERROR # 435; OPEN FILE FAILED FOR',A40)
      NERR=NERR+1

      END
C*********************************************************************
      SUBROUTINE RESIN (NERR)
C*********************************************************************

C  ROUTINE DIRECTS INPUT OF A RESTART FILE

C  NERR = ERROR COUNT STEPPED BY ONE ON ERROR (INPUT AND OUTPUT, INTEGER)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'
      INTEGER NERR
      INTEGER IU(6),I2,J2,N,I,J,K
      LOGICAL LSET

! bag8 - Under normal circumstances, set time information from restart
!        file, unless ZERO_RESTART flag is on.
      IF (ZERO_RESTART) THEN
        LSET=.FALSE.
      ELSE
        LSET=.TRUE.
      ENDIF

      IU(1)=0
      NERIN=NERR
      IF (LSET) THEN
        TIM=RTIMIN
      ELSE
        TIM=0.D0
      ENDIF

C  START INPUT LOOP

   51 IF (MYPRC.EQ.0) THEN
         IF (FORMIN) THEN
            READ (NFRIN,52) IU
   52       FORMAT(6I9)
         ELSE
            READ (NFRIN) IU
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREAD(6,IU)

      GOTO (1,2,3,4,5,6),IU(1)

      IF (LEVELC) WRITE (NFOUT,54) IU(1)
   54 FORMAT (' ERROR 436 - INVALID RESTART KEY =',I4)
      GO TO 13

C  INPUT RESTART CONTROL RECORD SET (L=IU(1)=1)

    1 IF (LSET) NSTEP=IU(2)
      I2=IU(3)
      J2=IU(4)
      N=NACTTIM-1
      IF (MYPRC.EQ.0) THEN
         WRITE (*,212) RFILNAM,NSTEP,TIM
  212    FORMAT(' READING RESTART FILE ',A40/
     & '      AT STEP',I6,' AND TIME',F12.3)
         IF (FORMIN) THEN
            READ (NFRIN,21) (RBUFR8(I),I=1,N)
   21       FORMAT(5G23.16)
            READ (NFRIN,21) (((BALANCE(I,J,K),I=1,I2),J=1,J2),K=1,8)
         ELSE
            READ (NFRIN) (RBUFR8(I),I=1,N)
            READ (NFRIN) (((BALANCE(I,J,K),I=1,I2),J=1,J2),K=1,8)
         ENDIF
      ENDIF
      IF (NUMPRC.GT.1) CALL SPREAD8(N,RBUFR8)
      IF (NUMPRC.GT.1) CALL SPREAD8(7*I2*J2,BALANCE)

      IF (LSET) THEN
      DELTIM=RBUFR8(1)
      DO 22 I=3,NACTTIM
   22 ACTTIM(I)=RBUFR8(I-1)
      IF (ACTTIM(3).LE.TIM) ACTTIM(3)=ACTTIM(3)+DTIMOUT
      ACTTIM(4)=TIM+DTIMRES
      ENDIF

      DO 23 I=1,I2
      DO 23 J=1,J2
      DO 23 K=1,3
   23 BALANCE(I,J,K)=0.D0

      GO TO 51

C  INPUT TERMINATION RECORD (L=2)

    2 NERR=NERIN
      CALL WAITALL()
      RETURN

C  INPUT RESTART ARRAYS (L=3)

    3 CALL RGEAIN(NERR)
      IF (NERR.GT.0) RETURN
      GO TO 51

C  INPUT WELL RESTART DATA (L=4)

    4 CONTINUE
      GO TO 51

C  INPUT LINEAR SOLVER DATA (L=6)

    6 CONTINUE
C      CALL LSORRI(IU,NERR)
      GO TO 51

C  INPUT PHYSICAL MODEL RESTART DATA (L=5)

    5 GO TO (61,62,63,64,65,66,67,60,69,60,60,60,73,60,60,76,77,78,79)
     &       ,IU(2)
   60 WRITE(*,*) ' ERROR - INVALID MODEL NUMBER IN RESTART HEADER'
      GO TO 13

   61 MODACT=1
C       CALL PRESTI(IU,NERR)
      GO TO 55
   62 MODACT=2
C       CALL IRESTI(IU,NERR)
      GO TO 55
   63 MODACT=3
C         CALL XRESTI(IU,NERR)
      GO TO 55
   64 MODACT=4
C         CALL CRESTI(IU,NERR)
      GO TO 55
   65 MODACT=5
C       CALL HRESTI(IU,NERR)
      GO TO 55
   66 MODACT=6
C       CALL GRESTI(IU,NERR)
      GO TO 55
   67 MODACT=7
C       CALL MRESTI(IU,NERR)
      GO TO 55
   69 MODACT=9
C      CALL SRESTI(IU,NERR)
      GO TO 55
   73 MODACT=13
      CALL TRESTI(IU,NERR)
      GO TO 55
   76 MODACT=16
C       CALL XRESTI(IU,NERR)
      GO TO 55
   77 MODACT=17
C      CALL TRESTI(IU,NERR)
      GO TO 55
   78 MODACT=18
C       CALL HRESTI(IU,NERR)
      GO TO 55
   79 MODACT=19
C       CALL HRESTI(IU,NERR)

   55 MODACT=0
      GO TO 51

   13 NERR=NERR+1
      RETURN

      END
C*********************************************************************
      SUBROUTINE RGEAOUT (NARR,NEXTRA,NERR)
C*********************************************************************

C  DIRECTS WRITING OF A GRID-ELEMENT ARRAY TO A RESTART FILE

C  NARR = GRID-ELEMENT ARRAY NUMBER OF ARRAY TO BE OUTPUT (INPUT, INTEGER)

C  NEXTRA = PRODUCT OF 4TH AND HIGHER DIMENSIONS TO BE OUTPUT (INPUT, INTEGER)
C           (SET TO 1 IF THERE IS NO FORTH DIMENSION)
C           NOTE THAT ONLY 1 VALUE OF THE 4TH DIMENSION IS OUTPUT FOR EACH
C           CALL TO RGEAOUT

C  NERR = ERROR COUNT STEPPED BY ONE ON ERROR (INPUT AND OUTPUT, INTEGER)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER NARG(2)
      CHARACTER*10 ANAM(2)
      EXTERNAL RSR4,RSR8,RSI4,RSL4

      MTM=CURRENT_MODEL+1
      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      CALL ARYTYPE(NARR,KIND,NDIM4,NMODIN,KERR)
      IF (KERR.GT.0) GO TO 13

      IF (MYPRC.EQ.0) THEN
         CALL GETANAM(NARR,ANAM(1),KERR)
         IF (KERR.GT.0) GO TO 13
         L=3
         IF (FORMOUT) THEN
            WRITE (NFROUT,11) L,NARR,L,L,L,L
   11       FORMAT(6I9)
            WRITE (NFROUT,12) ANAM(1)
   12       FORMAT(A10)
         ELSE
            WRITE (NFROUT) L,NARR,L,L,L,L
            WRITE (NFROUT) ANAM(1)
         ENDIF
      ENDIF

      IF (MYPRC.GT.0) CALL POLL(MSGTAG(MTM),NP)

      NARG(1)=1
      NARG(2)=NARR
      NARS=NARR
      NEORS=0
      I4UTIL=NEXTRA
      IF (I4UTIL.LT.1) I4UTIL=1
      IF (I4UTIL.GT.NDIM4) I4UTIL=NDIM4
      NERRP=NERR
      GOTO (1,2,3,3,5,5),KIND

    1 CALL CALLWORK(RSR4,NARG)
      IF (MYPRC.EQ.0) CALL RRR4(NERR)
      GO TO 10

    2 CALL CALLWORK(RSR8,NARG)
      IF (MYPRC.EQ.0) CALL RRR8(NERR)
      GO TO 10

    3 CALL CALLWORK(RSI4,NARG)
      IF (MYPRC.EQ.0) CALL RRI4(NERR)
      GO TO 10

    5 CALL CALLWORK(RSL4,NARG)
      IF (MYPRC.EQ.0) CALL RRI4(NERR)

   10 NERR=NERRP

      L=0
      IF (MYPRC.EQ.0) THEN
         IF (FORMOUT) THEN
            WRITE (NFROUT,15) L,L,L,L,L,L
   15       FORMAT(6I8)
         ELSE
            WRITE (NFROUT) L,L,L,L,L,L
         ENDIF
      ENDIF

      RETURN

   13 WRITE (NFOUT,14)
   14 FORMAT (/' ERROR # 434; INVALID ARRAY NUMBER IN RGEAOUT()')
      NERR=NERR+1

      END
C*********************************************************************
      SUBROUTINE RSR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                 KEYOUT,NBLK,ARY)
C*********************************************************************

C  WRITES A REAL*4 GRID-ELEMENT ARRAY TO A RESTART FILE

C  ARY = GRID-ELEMENT ARRAY TO BE OUTPUT (INPUT, REAL*4)

C  NOTE: 4TH DIMENSION IS PASSED IN I4UTIL IN /blkary/

C  NOTE: THE GRID-ELEMENT ARRAY NUMBER IS PASSED IN NARS IN /output/.
C        AFTER OUTPUT, NARS IS RESET TO -1 (PROCESSOR 0)

C  BUFFER:
C     H(1)=NBLK
C     H(2)=JG1
C     H(3)=JG2
C     H(4)=KG
C     H(5)=L   (4TH INDEX)
C     H(6)=N   (LAST H INDEX USED)

C*********************************************************************
      PARAMETER (NMX=3*480)
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'

      REAL*4   ARY(IDIM,JDIM,KDIM,*)
      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      LOGICAL  OUTA

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      JG1=0
      N=6
      RBUFR4(1)=NBLK
      RBUFR4(5)=I4UTIL

      DO 2 K=KL1,KL2
      JL1=JL1V(K)
      JL2=JL2V(K)
      KG=K+KOFF
      RBUFR4(4)=KG
      NK=N0MAP(NBLK)+KG*NYMAP(NBLK)+JOFF
      DO 2 J=JL1,JL2
      OUTA=.FALSE.

      IF (PRCMAP(NK+J).EQ.MYPRC) THEN
         JG2=J+JOFF
         IF (JG1.EQ.0) JG1=JG2
         DO 3 I=IL1,IL2
         N=N+1
    3    RBUFR4(N)=ARY(I,J,K,I4UTIL)
      ELSE
         OUTA=.TRUE.
      ENDIF

      IF (N+IDIM.GT.NMX) OUTA=.TRUE.
      IF (J.EQ.JL2) OUTA=.TRUE.
      IF (OUTA.AND.(N.GT.6)) THEN
         RBUFR4(2)=JG1
         RBUFR4(3)=JG2
         RBUFR4(6)=N

         IF (MYPRC.EQ.0) THEN

            NB=RBUFR4(1)+.1
            JA=RBUFR4(2)+.1
            JB=RBUFR4(3)+.1
            KK=RBUFR4(4)+.1
            ME=RBUFR4(5)+.1
            MM=RBUFR4(6)+.1
            NEORS=NEORS+MM-6
            IF (FORMOUT) THEN
               WRITE (NFROUT,4) NB,JA,JB,KK,ME,MM-6
    4          FORMAT(6I8)
               WRITE (NFROUT,5) (RBUFR4(NN),NN=7,MM)
    5          FORMAT(8G15.8)
            ELSE
               WRITE (NFROUT) NB,JA,JB,KK,ME,MM-6
               WRITE (NFROUT) (RBUFR4(NN),NN=7,MM)
            ENDIF
         ELSE

      CALL RSSR4(N,RBUFR4,NERR)

         ENDIF
         N=6
         JG1=0
      ENDIF

    2 CONTINUE

      END
C*********************************************************************
      SUBROUTINE RSR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                 KEYOUT,NBLK,ARY)
C*********************************************************************

C  WRITES A REAL*8 GRID-ELEMENT ARRAY TO A RESTART FILE

C  ARY = GRID-ELEMENT ARRAY TO BE OUTPUT (INPUT, REAL*8)

C*********************************************************************
      PARAMETER (NMX=3*480)
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'

      REAL*8   ARY(IDIM,JDIM,KDIM,*)
      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      LOGICAL  OUTA

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      JG1=0
      N=6
      RBUFR8(1)=NBLK
      RBUFR8(5)=I4UTIL

      DO 2 K=KL1,KL2
      JL1=JL1V(K)
      JL2=JL2V(K)
      KG=K+KOFF
      RBUFR8(4)=KG
      NK=N0MAP(NBLK)+KG*NYMAP(NBLK)+JOFF
      DO 2 J=JL1,JL2
      OUTA=.FALSE.

      IF (PRCMAP(NK+J).EQ.MYPRC) THEN
         JG2=J+JOFF
         IF (JG1.EQ.0) JG1=JG2
         DO 3 I=IL1,IL2
         N=N+1
    3    RBUFR8(N)=ARY(I,J,K,I4UTIL)
      ELSE
         OUTA=.TRUE.
      ENDIF

      IF (N+IDIM.GT.NMX) OUTA=.TRUE.
      IF (J.EQ.JL2) OUTA=.TRUE.
      IF (OUTA.AND.(N.GT.6)) THEN
         RBUFR8(2)=JG1
         RBUFR8(3)=JG2
         RBUFR8(6)=N

         IF (MYPRC.EQ.0) THEN

            NB=RBUFR8(1)+.1
            JA=RBUFR8(2)+.1
            JB=RBUFR8(3)+.1
            KK=RBUFR8(4)+.1
            ME=RBUFR8(5)+.1
            MM=RBUFR8(6)+.1
            NEORS=NEORS+MM-6
            IF (FORMOUT) THEN
               WRITE (NFROUT,4) NB,JA,JB,KK,ME,MM-6
    4          FORMAT(6I8)
               WRITE (NFROUT,5) (RBUFR8(NN),NN=7,MM)
    5          FORMAT(5G23.16)
            ELSE
               WRITE (NFROUT) NB,JA,JB,KK,ME,MM-6
               WRITE (NFROUT) (RBUFR8(NN),NN=7,MM)
            ENDIF

         ELSE

      CALL RSSR8(N,RBUFR8,NERR)

         ENDIF
         N=6
         JG1=0
      ENDIF

    2 CONTINUE

      END
C*********************************************************************
      SUBROUTINE RSI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                 KEYOUT,NBLK,IARY)
C*********************************************************************

C  WRITES A INTEGER GRID-ELEMENT ARRAY TO A RESTART FILE

C  IARY = GRID-ELEMENT ARRAY TO BE OUTPUT (INPUT, INTEGER)

C*********************************************************************
      PARAMETER (NMX=3*480)
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'

      INTEGER  IARY(IDIM,JDIM,KDIM,*)
      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      LOGICAL  OUTA

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      JG1=0
      N=6
      RBUFI4(1)=NBLK
      RBUFI4(5)=I4UTIL

      DO 2 K=KL1,KL2
      JL1=JL1V(K)
      JL2=JL2V(K)
      KG=K+KOFF
      RBUFI4(4)=KG
      NK=N0MAP(NBLK)+KG*NYMAP(NBLK)+JOFF
      DO 2 J=JL1,JL2
      OUTA=.FALSE.

      IF (PRCMAP(NK+J).EQ.MYPRC) THEN
         JG2=J+JOFF
         IF (JG1.EQ.0) JG1=JG2
         DO 3 I=IL1,IL2
         N=N+1
    3    RBUFI4(N)=IARY(I,J,K,I4UTIL)
      ELSE
         OUTA=.TRUE.
      ENDIF

      IF (N+IDIM.GT.NMX) OUTA=.TRUE.
      IF (J.EQ.JL2) OUTA=.TRUE.
      IF (OUTA.AND.(N.GT.6)) THEN
         RBUFI4(2)=JG1
         RBUFI4(3)=JG2
         RBUFI4(6)=N

         IF (MYPRC.EQ.0) THEN

            MM=RBUFI4(6)
            NEORS=NEORS+MM-6
            IF (FORMOUT) THEN
               WRITE (NFROUT,4) (RBUFI4(NN),NN=1,5),MM-6
    4          FORMAT(6I8)
               WRITE (NFROUT,5) (RBUFI4(NN),NN=7,MM)
    5          FORMAT(12I9)
            ELSE
               WRITE (NFROUT) (RBUFI4(NN),NN=1,5),MM-6
               WRITE (NFROUT) (RBUFI4(NN),NN=7,MM)
            ENDIF

         ELSE

      CALL RSSI4(N,RBUFI4,NERR)

         ENDIF
         N=6
         JG1=0
      ENDIF

    2 CONTINUE

      END
C*********************************************************************
      SUBROUTINE RSL4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                 KEYOUT,NBLK,LARY)
C*********************************************************************

C  WRITES A LOGICAL GRID-ELEMENT ARRAY TO A RESTART FILE

C  LARY = GRID-ELEMENT ARRAY TO BE OUTPUT (INPUT, LOGICAL)

C*********************************************************************
      PARAMETER (NMX=3*480)
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'blkary.h'
      INCLUDE 'restc.h'

      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      LOGICAL  OUTA,LARY(IDIM,JDIM,KDIM,*)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      JG1=0
      N=6
      RBUFI4(1)=NBLK
      RBUFI4(5)=I4UTIL

      DO 2 K=KL1,KL2
      JL1=JL1V(K)
      JL2=JL2V(K)
      KG=K+KOFF
      RBUFI4(4)=KG
      NK=N0MAP(NBLK)+KG*NYMAP(NBLK)+JOFF
      DO 2 J=JL1,JL2
      OUTA=.FALSE.

      IF (PRCMAP(NK+J).EQ.MYPRC) THEN
         JG2=J+JOFF
         IF (JG1.EQ.0) JG1=JG2
         DO 3 I=IL1,IL2
         N=N+1
         IF (LARY(I,J,K,I4UTIL)) THEN
            RBUFI4(N)=1
         ELSE
            RBUFI4(N)=0
         ENDIF
    3    CONTINUE
      ELSE
         OUTA=.TRUE.
      ENDIF

      IF (N+IDIM.GT.NMX) OUTA=.TRUE.
      IF (J.EQ.JL2) OUTA=.TRUE.
      IF (OUTA.AND.(N.GT.6)) THEN
         RBUFI4(2)=JG1
         RBUFI4(3)=JG2
         RBUFI4(6)=N

         IF (MYPRC.EQ.0) THEN

            MM=RBUFI4(6)
            NEORS=NEORS+MM-6
            IF (FORMOUT) THEN
               WRITE (NFROUT,4) (RBUFI4(NN),NN=1,5),MM-6
    4          FORMAT(6I8)
               WRITE (NFROUT,5) (RBUFI4(NN),NN=7,MM)
    5          FORMAT(55I2)
            ELSE
               WRITE (NFROUT) (RBUFI4(NN),NN=1,5),MM-6
               WRITE (NFROUT) (RBUFI4(NN),NN=7,MM)
            ENDIF

         ELSE

      CALL RSSI4(N,RBUFI4,NERR)

         ENDIF
         N=6
         JG1=0
      ENDIF

    2 CONTINUE

      END

C*********************************************************************
      SUBROUTINE RGEAIN (NERR)
C*********************************************************************

C  DIRECTS INPUT OF A GRID-ELEMENT ARRAY FROM A RESTART FILE

C*********************************************************************
      PARAMETER (NMX=3*480)
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER NARG(2)
      CHARACTER*10 ANAM(2)

      EXTERNAL RSINR4,RSINR8,RSINI4,RSINL4

      IF (MYPRC.EQ.0) THEN
         IF (FORMIN) THEN
            READ (NFRIN,1) ANAM(1)
    1       FORMAT(A10)
         ELSE
            READ (NFRIN) ANAM(1)
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREADC(10,ANAM(1))

      CALL ARYDAT(ANAM(1),KINDIN,ND4IN,NAIN,IERR)
      CALL ARYTYPE(NAIN,KDUM1,KDUM2,NMODIN,IERR)
      NARG(1)=1
      NARG(2)=NAIN

      MTM=CURRENT_MODEL+1
      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      NERIN=NERR

    4 IF (MYPRC.EQ.0) THEN
         IF (FORMIN) THEN
            READ (NFRIN,2) NBIN,JAIN,JBIN,KIN,MEIN,MMIN
    2       FORMAT(6I8)
         ELSE
            READ (NFRIN) NBIN,JAIN,JBIN,KIN,MEIN,MMIN
         ENDIF
      ENDIF

      GOTO (11,12,13,13,15,15),KINDIN

C  -----------------------  REAL*4

   11 IF (MYPRC.EQ.0) THEN
         RBUFR4(1)=NBIN
         RBUFR4(2)=JAIN
         RBUFR4(3)=JBIN
         RBUFR4(4)=KIN
         RBUFR4(5)=MEIN
         RBUFR4(6)=MMIN
         MM=MMIN+6
         IF (MM.GT.6) THEN
            IF (FORMIN) THEN
               READ (NFRIN,21) (RBUFR4(NN),NN=7,MM)
   21          FORMAT(8G15.8)
            ELSE
               READ (NFRIN) (RBUFR4(NN),NN=7,MM)
            ENDIF
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREAD (1,MM)
      IF (NUMPRC.GT.1) CALL SPREAD4 (MM,RBUFR4)

      NBIN=RBUFR4(1)+.1

      IF (NBIN.EQ.0) RETURN

      JAIN=RBUFR4(2)+.1
      JBIN=RBUFR4(3)+.1
      KIN=RBUFR4(4)+.1
      MEIN=RBUFR4(5)+.1
      MMIN=RBUFR4(6)+.1

      MODACT=NMODIN
      CALL CALLWORK(RSINR4,NARG)
      MODACT=0
      GO TO 5

C  -----------------------  REAL*8

   12 IF (MYPRC.EQ.0) THEN
         RBUFR8(1)=NBIN
         RBUFR8(2)=JAIN
         RBUFR8(3)=JBIN
         RBUFR8(4)=KIN
         RBUFR8(5)=MEIN
         RBUFR8(6)=MMIN
         MM=MMIN+6
         IF (MM.GT.6) THEN
            IF (FORMIN) THEN
               READ (NFRIN,22) (RBUFR8(NN),NN=7,MM)
   22          FORMAT(5G23.16)
            ELSE
               READ (NFRIN) (RBUFR8(NN),NN=7,MM)
            ENDIF
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREAD (1,MM)
      IF (NUMPRC.GT.1) CALL SPREAD8 (MM,RBUFR8)

      NBIN=RBUFR8(1)+.1

      IF (NBIN.EQ.0) RETURN

      JAIN=RBUFR8(2)+.1
      JBIN=RBUFR8(3)+.1
      KIN=RBUFR8(4)+.1
      MEIN=RBUFR8(5)+.1
      MMIN=RBUFR8(6)+.1

      MODACT=NMODIN
      CALL CALLWORK(RSINR8,NARG)
      MODACT=0
      GO TO 5

C  -----------------------  INTEGER*4

   13 IF (MYPRC.EQ.0) THEN
         RBUFI4(1)=NBIN
         RBUFI4(2)=JAIN
         RBUFI4(3)=JBIN
         RBUFI4(4)=KIN
         RBUFI4(5)=MEIN
         RBUFI4(6)=MMIN
         MM=MMIN+6
         IF (MM.GT.6) THEN
            IF (FORMIN) THEN
               READ (NFRIN,23) (RBUFI4(NN),NN=7,MM)
   23          FORMAT(12I9)
            ELSE
               READ (NFRIN) (RBUFI4(NN),NN=7,MM)
            ENDIF
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREAD (1,MM)
      IF (NUMPRC.GT.1) CALL SPREAD (MM,RBUFI4)

      NBIN=RBUFI4(1)

      IF (NBIN.EQ.0) RETURN

      JAIN=RBUFI4(2)
      JBIN=RBUFI4(3)
      KIN=RBUFI4(4)
      MEIN=RBUFI4(5)
      MMIN=RBUFI4(6)

      MODACT=NMODIN
      CALL CALLWORK(RSINI4,NARG)
      MODACT=0
      GO TO 5

C  -----------------------  LOGICAL*4

   15 IF (MYPRC.EQ.0) THEN
         RBUFI4(1)=NBIN
         RBUFI4(2)=JAIN
         RBUFI4(3)=JBIN
         RBUFI4(4)=KIN
         RBUFI4(5)=MEIN
         RBUFI4(6)=MMIN
         MM=MMIN+6
         IF (MM.GT.6) THEN
            IF (FORMIN) THEN
               READ (NFRIN,25) (RBUFI4(NN),NN=7,MM)
   25          FORMAT(55I2)
            ELSE
               READ (NFRIN) (RBUFI4(NN),NN=7,MM)
            ENDIF
         ENDIF
      ENDIF

      IF (NUMPRC.GT.1) CALL SPREAD (1,MM)
      IF (NUMPRC.GT.1) CALL SPREAD (MM,RBUFI4)

      NBIN=RBUFI4(1)

      IF (NBIN.EQ.0) RETURN

      JAIN=RBUFI4(2)
      JBIN=RBUFI4(3)
      KIN=RBUFI4(4)
      MEIN=RBUFI4(5)
      MMIN=RBUFI4(6)

      MODACT=NMODIN
      CALL CALLWORK(RSINL4,NARG)
      MODACT=0

    5 NERR=NERIN
      GO TO 4

      END

C*********************************************************************
      SUBROUTINE RSINR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,ARY)
C*********************************************************************

C  INSERTS REAL*4 ARRAY VALUES DURING RESTART INPUT

C  ARY = GRID-ELEMENT ARRAY TO BE INPUT (INPUT, REAL*4)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)

      REAL*4 ARY(IDIM,JDIM,KDIM,*)

      IF (NBIN.NE.NBLK) RETURN

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      K=KIN-KOFF
      IF (K.LT.KL1.OR.K.GT.KL2) RETURN

      N=6
      J1=JAIN-JOFF
      JL1=JL1V(K)
      J2=JBIN-JOFF
      JL2=JL2V(K)
      NI=IL2-IL1+1

      DO 1 J=J1,J2
      IF (J.GE.JL1.AND.J.LE.JL2) THEN
         DO 2 I=IL1,IL2
         N=N+1
    2    ARY(I,J,K,MEIN)=RBUFR4(N)
      ELSE
         N=N+NI
      ENDIF
    1 CONTINUE

      END

C*********************************************************************
      SUBROUTINE RSINR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,ARY)
C*********************************************************************

C  INSERTS REAL*8 ARRAY VALUES DURING RESTART INPUT

C  ARY = GRID-ELEMENT ARRAY TO BE INPUT (INPUT, REAL*8)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)

      REAL*8 ARY(IDIM,JDIM,KDIM,*)

      IF (NBIN.NE.NBLK) RETURN

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      K=KIN-KOFF
      IF (K.LT.KL1.OR.K.GT.KL2) RETURN

      N=6
      J1=JAIN-JOFF
      JL1=JL1V(K)
      J2=JBIN-JOFF
      JL2=JL2V(K)
      NI=IL2-IL1+1

      DO 1 J=J1,J2
      IF (J.GE.JL1.AND.J.LE.JL2) THEN
         DO 2 I=IL1,IL2
         N=N+1
    2    ARY(I,J,K,MEIN)=RBUFR8(N)
      ELSE
         N=N+NI
      ENDIF
    1 CONTINUE

      END

C*********************************************************************
      SUBROUTINE RSINI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,IARY)
C*********************************************************************

C  INSERTS INTEGER ARRAY VALUES DURING RESTART INPUT

C  IARY = GRID-ELEMENT ARRAY TO BE INPUT (INPUT, INTEGER)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)

      INTEGER IARY(IDIM,JDIM,KDIM,*)

      IF (NBIN.NE.NBLK) RETURN

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      K=KIN-KOFF
      IF (K.LT.KL1.OR.K.GT.KL2) RETURN

      N=6
      J1=JAIN-JOFF
      JL1=JL1V(K)
      J2=JBIN-JOFF
      JL2=JL2V(K)
      NI=IL2-IL1+1

      DO 1 J=J1,J2
      IF (J.GE.JL1.AND.J.LE.JL2) THEN
         DO 2 I=IL1,IL2
         N=N+1
    2    IARY(I,J,K,MEIN)=RBUFI4(N)
      ELSE
         N=N+NI
      ENDIF
    1 CONTINUE

      END

C*********************************************************************
      SUBROUTINE RSINL4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,IARY)
C*********************************************************************

C  INSERTS LOGICAL ARRAY VALUES DURING RESTART INPUT

C  IARY = GRID-ELEMENT ARRAY TO BE INPUT (INPUT, LOGICAL)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'output.h'
      INCLUDE 'restc.h'

      INTEGER  JL1V(KDIM),JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)

      LOGICAL IARY(IDIM,JDIM,KDIM,*)

      IF (NBIN.NE.NBLK) RETURN

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      K=KIN-KOFF
      IF (K.LT.KL1.OR.K.GT.KL2) RETURN

      N=6
      J1=JAIN-JOFF
      JL1=JL1V(K)
      J2=JBIN-JOFF
      JL2=JL2V(K)
      NI=IL2-IL1+1

      DO 1 J=J1,J2
      IF (J.GE.JL1.AND.J.LE.JL2) THEN
         DO 2 I=IL1,IL2
         N=N+1
         IF (RBUFI4(N).EQ.0) THEN
            IARY(I,J,K,MEIN)=.FALSE.
         ELSE
            IARY(I,J,K,MEIN)=.TRUE.
         ENDIF
    2    CONTINUE
      ELSE
         N=N+NI
      ENDIF
    1 CONTINUE

      END
