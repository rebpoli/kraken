C  MANY.F - MULTIPROCESSOR ROUTINES

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE SETPRCS (NERR)
C  SUBROUTINE KILLPRC (NERR)
C  SUBROUTINE SNDABUF (PKEY)
C  SUBROUTINE COMMI   (NERR)
C  SUBROUTINE EXMLST  ()
C  SUBROUTINE R4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE R8UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE I4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE I2UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE L4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE L2UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C  SUBROUTINE NOUPDATE (NBLK,IV1,IV2)
C  SUBROUTINE R4TOBUF (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE BUFTOR4 (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE R8TOBUF (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE BUFTOR8 (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE I4TOBUF (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE BUFTOI4 (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE I2TOBUF (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE BUFTOI2 (IDIM,JDIM,KDIM,AR)
C  SUBROUTINE GETMP (I1G,I2G,ITS,IOFF,J1G,J2G,JTS,JOFF,KG,KOFF,
C                    IEX,KEYOUT,IDIM,JDIM,KDIM,NBLK,KARY,KIND,ARRY)
C  SUBROUTINE R4GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY)
C  SUBROUTINE R8GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY)
C  SUBROUTINE I2GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY)
C  SUBROUTINE I4GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY)

C  CODE HISTORY:

C  JOHN WHEELER    12/14/95     ALPHA CODE
C  JOHN WHEELER     6/22/96     DEBUG
C  JOHN WHEELER     7/12/98     IMPROVE UPDATE EFFICIENCY
C  JOHN WHEELER     4/20/99     ADD GENERIC MULTI-MODEL CAPABILITY

C  NOTES:

C     1)  The routines in this file are required when MANY is invoked
C         in the .siz file.

C     2)  ERROR NUMBERS 201 TO 240 ARE RESERVED FOR PARALLEL ROUTINES

C*********************************************************************
      SUBROUTINE SETPRCS (NERR)
C*********************************************************************

C  Routine sets multiprocessor parameters including number of processors,
C  processor number, and process id (if appropriate).

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER NERR,IERR
      CHARACTER*14 PARROU

      IERR=0
      PARROU='MPI_INIT'
! dealii - check if MPI already initialized for using ipars as a library
      CALL MPI_INITIALIZED(MPI_EXTERNAL,IERR)
      IF (.NOT.MPI_EXTERNAL) THEN
        CALL MPI_INIT (IERR)
        IF (IERR.GT.0) GO TO 13
      ENDIF
      PARROU='MPI_COMM_RANK'
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPRC,IERR)
      IF (IERR.GT.0) GO TO 13
      PARROU='MPI_COMM_SIZE'
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NUMPRC,IERR)
      IF (IERR.GT.0) GO TO 13

! bag8 - error if size parameters aren't big enough
      IF (NUMPRC.GT.256) THEN
        IF (MYPRC.EQ.0)
     &    WRITE(*,'(2(A,I))')'Error in SETPRCS: NUMPRC=',NUMPRC,
     &      'greater than MXPROC=',256
        CALL KILLPRC(1)
      ENDIF

      RETURN

   13 NERR=NERR+1
      WRITE (*,14) PARROU
   14 FORMAT (/' ERROR # 202; PARALLEL ROUTINE ',A14,' FAILED')

      END

C*********************************************************************
      SUBROUTINE KILL_IPARS(MSG)
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      CHARACTER*(*) :: MSG
      IF (MYPRC.EQ.0) WRITE(*,'(2A)') 'KILL_IPARS: ',MSG
      CALL KILLPRC(1)
      END

C*********************************************************************
      SUBROUTINE KILLPRC (NERR)
C*********************************************************************

C  Routine terminates a multiprocessor simulation.

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER NERR,IERR

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE KILLPRC'

      CALL VISQUIT(NERR)
      IF (GEADBG.GT.0) CALL CLOSE_GEA()
      IF (NERR.EQ.0) CALL WAITALL()

      IERR=0
      IF (.NOT.MPI_EXTERNAL) CALL MPI_FINALIZE(IERR)

      IF (NERR.EQ.0) THEN
C         IF (DEALII) RETURN
         STOP 0
      ELSE
         STOP 13
      ENDIF

      END
C*********************************************************************
      SUBROUTINE SNDABUF (PKEY)
C*********************************************************************

C  Sends a buffer of character data from processor 0 to all other
C  processors.  A blocking broadcast of the A() buffer in /SCRAT1/ is
C  used.  LAST in /SCRAT1/ is set to the number of characters received.
C  This routine can be called only by the framework.

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'mpif.h'
      LOGICAL PKEY

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE SNDABUF'

      CALL TIMON(2)

      CALL MPI_BCAST(LAST,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL WAITALL()
      CALL MPI_BCAST(A,LAST,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

      CALL TIMOFF(2)

      IF (PKEY) THEN
         IF (MYPRC.EQ.0) WRITE (*,*)
     &     'DATA BROADCAST FROM PROCESSOR 0,',LAST,' BYTES'
         IF (MYPRC.EQ.1) WRITE (*,*)
     &     'BROADCAST DATA RECEIVED AT PROCESSOR 1,',LAST,' BYTES'
      ENDIF

      END
C*********************************************************************
      SUBROUTINE COMMI (NERR)
C*********************************************************************

C  Sets up message system for updating the communication layer around
C  the elements assigned to a processor.  This routine should not be
C  called for dynamic load balancing.  Difference templates currently
C  supported are as follows:

C  n = 1 ==> adjacent elements only (7 point template)
C    = 2 ==> cube (27 point template) (also cornerless cube)
C    = 3 ==> "higher order" (27 point template)

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************

      INCLUDE 'control.h'
      INCLUDE 'layout.h'

      INTEGER JS1(4),KS1(4),JS2(8),KS2(8),JS3(12),KS3(12)
      DATA JS1/1,-1,0,0/,KS1/0,0,1,-1/
      DATA JS2/-1,-1,-1,0,0,1,1,1/,KS2/-1,0,1,-1,1,-1,0,1/
      DATA JS3/-1,-1,-1,0,0,1,1,1,-2,0,0,2/,
     &     KS3/-1,0,1,-1,1,-1,0,1,0,2,-2,0/
      INTEGER NT1,NT2,NT3

C  INITIALIZE

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,*)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE COMMI'

      LSTIND=0
      MXTP=4
      IF (MXTP.LT.1) RETURN
      DO 100 I=1,MXTP
      DO 100 J=1,10
  100 MSGS(I,J)=0
      IF (JLAY.LT.1.OR.KLAY.LT.1) RETURN

C  7 POINT TEMPLATE

      NT1=1
      NUMTMP=1
      DO 1 M=1,NUMBLK
      CALL BLKDIM (M,NX,NY,NZ,KERR)
      N0=N0MAP(M)
      NYM=NYMAP(M)
      MG=0

C  BUILD LIST OF PROCESSOR TARGETS FOR 7 POINT TEMPLATE

      DO 2 J=1,NY
      DO 2 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 3 MM=1,4
         JT=J+JS1(MM)
         KT=K+KS1(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MR=PRCMAP(N0+NYM*KT+JT)
            IF (MR.NE.MYPRC.AND.MR.GE.0) THEN
               DO 4 L=1,MG
               IF (MR.EQ.MSGTRG(L,NT1,M)) GO TO 3
    4          CONTINUE
               IF (MG.EQ.10) GO TO 213
               MG=MG+1
               MSGTRG(MG,NT1,M)=MR
            ENDIF
         ENDIF
    3    CONTINUE
      ENDIF
    2 CONTINUE
      MSGS(NT1,M)=MG

C  BUILD MESSAGE LISTS FOR 7 POINT TEMPLATE

      DO 5 MN=1,MG
      MR=MSGTRG(MN,NT1,M)
      MSGSNDF(MN,NT1,M)=LSTIND+1
      MSGSNDL(MN,NT1,M)=LSTIND
      DO 5 J=1,NY
      DO 5 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 6 MM=1,4
         JT=J+JS1(MM)
         KT=K+KS1(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MF=PRCMAP(N0+NYM*KT+JT)
            IF (MF.EQ.MR.AND.MF.GE.0) THEN
               MJK1=MSGSNDF(MN,NT1,M)
               MJK2=MSGSNDL(MN,NT1,M)
               DO 7 MJK=MJK1,MJK2
               IF (MSGIND(1,MJK).EQ.J.AND.MSGIND(2,MJK).EQ.K) GO TO 6
    7          CONTINUE
               IF (LSTIND.GE.2000) GO TO 113
               LSTIND=LSTIND+1
               MSGSNDL(MN,NT1,M)=LSTIND
               MSGIND(1,LSTIND)=J
               MSGIND(2,LSTIND)=K
            ENDIF
         ENDIF
    6    CONTINUE
      ENDIF
    5 CONTINUE

    1 CONTINUE

C  CUBE TEMPLATE

      NT2=2         ! bag8: NT2 defined to allow compilation with check bounds
      IF (MXTP.LT.NT2) GO TO 80
      NUMTMP=NUMTMP+1
      DO 11 M=1,NUMBLK
      CALL BLKDIM (M,NX,NY,NZ,KERR)
      N0=N0MAP(M)
      NYM=NYMAP(M)
      MG=0

C  BUILD LIST OF PROCESSOR TARGETS FOR CUBE TEMPLATE

      DO 12 J=1,NY
      DO 12 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 13 MM=1,8
         JT=J+JS2(MM)
         KT=K+KS2(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MR=PRCMAP(N0+NYM*KT+JT)
            IF (MR.NE.MYPRC.AND.MR.GE.0) THEN
               DO 14 L=1,MG
               IF (MR.EQ.MSGTRG(L,NT2,M)) GO TO 13
   14          CONTINUE
               IF (MG.EQ.10) GO TO 213
               MG=MG+1
               MSGTRG(MG,NT2,M)=MR
            ENDIF
         ENDIF
   13    CONTINUE
      ENDIF
   12 CONTINUE
      MSGS(NT2,M)=MG

C  BUILD MESSAGE LISTS FOR CUBE TEMPLATE

      DO 15 MN=1,MG
      MR=MSGTRG(MN,NT2,M)
      MSGSNDF(MN,NT2,M)=LSTIND+1
      MSGSNDL(MN,NT2,M)=LSTIND
      DO 15 J=1,NY
      DO 15 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 16 MM=1,8
         JT=J+JS2(MM)
         KT=K+KS2(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MF=PRCMAP(N0+NYM*KT+JT)
            IF (MF.EQ.MR.AND.MF.GE.0) THEN
               MJK1=MSGSNDF(MN,NT2,M)
               MJK2=MSGSNDL(MN,NT2,M)
               DO 17 MJK=MJK1,MJK2
               IF (MSGIND(1,MJK).EQ.J.AND.MSGIND(2,MJK).EQ.K) GO TO 16
   17          CONTINUE
               IF (LSTIND.GE.2000) GO TO 113
               LSTIND=LSTIND+1
               MSGSNDL(MN,NT2,M)=LSTIND
               MSGIND(1,LSTIND)=J
               MSGIND(2,LSTIND)=K
            ENDIF
         ENDIF
   16    CONTINUE
      ENDIF
   15 CONTINUE

   11 CONTINUE

C  HIGHER ORDER TEMPLATE

      NT3=3
      IF (MXTP.LT.NT3.OR.JLAY.LT.2.OR.KLAY.LT.2) GO TO 80
      NUMTMP=NUMTMP+1
      DO 21 M=1,NUMBLK
      CALL BLKDIM (M,NX,NY,NZ,KERR)
      N0=N0MAP(M)
      NYM=NYMAP(M)
      MG=0

C  BUILD LIST OF PROCESSOR TARGETS FOR HIGHER ORDER TEMPLATE

      DO 22 J=1,NY
      DO 22 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 23 MM=1,12
         JT=J+JS3(MM)
         KT=K+KS3(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MR=PRCMAP(N0+NYM*KT+JT)
            IF (MR.NE.MYPRC.AND.MR.GE.0) THEN
               DO 24 L=1,MG
               IF (MR.EQ.MSGTRG(L,NT3,M)) GO TO 23
   24          CONTINUE
               IF (MG.EQ.10) GO TO 213
               MG=MG+1
               MSGTRG(MG,NT3,M)=MR
            ENDIF
         ENDIF
   23    CONTINUE
      ENDIF
   22 CONTINUE
      MSGS(NT3,M)=MG

C  BUILD MESSAGE LISTS FOR HIGHER ORDER TEMPLATE

      DO 25 MN=1,MG
      MR=MSGTRG(MN,NT3,M)
      MSGSNDF(MN,NT3,M)=LSTIND+1
      MSGSNDL(MN,NT3,M)=LSTIND
      DO 25 J=1,NY
      DO 25 K=1,NZ
      MS=PRCMAP(N0+NYM*K+J)
      IF (MS.EQ.MYPRC) THEN
         DO 26 MM=1,12
         JT=J+JS3(MM)
         KT=K+KS3(MM)
         IF (JT.GT.0.AND.JT.LE.NY.AND.KT.GT.0.AND.KT.LE.NZ) THEN
            MF=PRCMAP(N0+NYM*KT+JT)
            IF (MF.EQ.MR.AND.MF.GE.0) THEN
               MJK1=MSGSNDF(MN,NT3,M)
               MJK2=MSGSNDL(MN,NT3,M)
               DO 27 MJK=MJK1,MJK2
               IF (MSGIND(1,MJK).EQ.J.AND.MSGIND(2,MJK).EQ.K) GO TO 26
   27          CONTINUE
               IF (LSTIND.GE.2000) GO TO 113
               LSTIND=LSTIND+1
               MSGSNDL(MN,NT3,M)=LSTIND
               MSGIND(1,LSTIND)=J
               MSGIND(2,LSTIND)=K
            ENDIF
         ENDIF
   26    CONTINUE
      ENDIF
   25 CONTINUE

   21 CONTINUE

CBW
!bw C CONSTRUCT PRCMAPN() FOR NODAL BASED PROCESSOR ASSIGNMENT
!bw       DO M=1,NUMBLK
!bw          CALL BLKDIM(M,NX,NY,NZ,KERR)
!bw          N0=N0MAP(M)
!bw          N0N=N0MAPN(M)
!bw          NYM=NYMAP(M)
!bw          NYMN=NYM+1
!bw          DO K=1,NZ
!bw             DO J=1,NY
!bw                IF (PRCMAP(N0+K*NYM+J).GE.0) THEN
!bw                   PRCMAPN(N0N+K*NYMN+J)=PRCMAP(N0+K*NYM+J)
!bw                   PRCMAPN(N0N+K*NYMN+J+1)=PRCMAP(N0+K*NYM+J)
!bw                   PRCMAPN(N0N+(K+1)*NYMN+J)=PRCMAP(N0+K*NYM+J)
!bw                   PRCMAPN(N0N+(K+1)*NYMN+J+1)=PRCMAP(N0+K*NYM+J)
!bw                ENDIF
!bw             ENDDO
!bw          ENDDO
!bw       ENDDO

C  NODAL BASED TEMPLATE
      IF (MXTP.LT.4.OR.JLAY.LT.2.OR.KLAY.LT.2) GO TO 80
      NUMTMP=NUMTMP+1
      DO M=1,NUMBLK
         CALL BLKDIM (M,NX,NY,NZ,KERR)
         N0N=N0MAPN(M)
         NYMN=NYMAP(M)+1
         MG=0

C  BUILD LIST OF PROCESSOR TARGETS FOR NODAL BASED TEMPLATE

         DO J=1,NY+1
            DO K=1,NZ+1
            MS=PRCMAPN(N0N+NYMN*K+J)
            IF (MS.EQ.MYPRC) THEN
               DO MM=1,8
                  JT=J+JS2(MM)
                  KT=K+KS2(MM)
                  IF (JT.GT.0.AND.JT.LE.(NY+1).AND.KT.GT.0.AND.
     &               KT.LE.(NZ+1)) THEN
                     MR=PRCMAPN(N0N+NYMN*KT+JT)
                     IF (MR.NE.MYPRC.AND.MR.GE.0) THEN
                        DO L=1,MG
                           IF (MR.EQ.MSGTRG(L,4,M)) GOTO 33
                        ENDDO
                        IF (MG.EQ.10) GO TO 213
                        MG=MG+1
                        MSGTRG(MG,4,M)=MR
                     ENDIF
                  ENDIF
  33              CONTINUE
               ENDDO
            ENDIF
            ENDDO
         ENDDO
         MSGS(4,M)=MG

C  BUILD MESSAGE LISTS FOR NODAL BASED TEMPLATE

         DO MN=1,MG
            MR=MSGTRG(MN,4,M)
            MSGSNDF(MN,4,M)=LSTIND+1
            MSGSNDL(MN,4,M)=LSTIND
            DO J=1,NY+1
               DO K=1,NZ+1
                  MS=PRCMAPN(N0N+NYMN*K+J)
                  IF (MS.EQ.MYPRC) THEN
                     DO MM=1,8
                        JT=J+JS2(MM)
                        KT=K+KS2(MM)
                        IF (JT.GT.0.AND.JT.LE.(NY+1).AND.
     &                      KT.GT.0.AND.KT.LE.(NZ+1)) THEN
                           MF=PRCMAPN(N0N+NYMN*KT+JT)
                           IF (MF.EQ.MR.AND.MF.GE.0) THEN
                              MJK1=MSGSNDF(MN,4,M)
                              MJK2=MSGSNDL(MN,4,M)
                              DO MJK=MJK1,MJK2
                                 IF (MSGIND(1,MJK).EQ.J.AND.
     &                               MSGIND(2,MJK).EQ.K) GOTO 36
                              ENDDO
                              IF (LSTIND.GE.2000) GO TO 113
                              LSTIND=LSTIND+1
                              MSGSNDL(MN,4,M)=LSTIND
                              MSGIND(1,LSTIND)=J
                              MSGIND(2,LSTIND)=K
                           ENDIF
                        ENDIF
 36                     CONTINUE
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
CBW

C  EXCHANGE MESSAGE LISTS

   80 CALL EXMLST ()

      RETURN

  113 NERR=NERR+1
      WRITE (*,114) MYPRC
  114 FORMAT (/' ERROR # 201; TOO MANY MESSAGE ELEMENTS, PROC',I5)
      RETURN

  213 NERR=NERR+1
      WRITE (*,214) MYPRC
  214 FORMAT (/' ERROR # 201; TOO MANY MESSAGE TARGETS, PROC',I5)

      END

! bag8
C*********************************************************************
      SUBROUTINE CHECK_COMMI()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'layout.h'
      IF (.NOT.COMMI_DONE) THEN
        CALL KILL_IPARS('Error: doing UPDATE before COMMI')
      ENDIF
      END

C*********************************************************************
      SUBROUTINE EXMLST ()
C*********************************************************************

C  Exchanges message lists for updating the communication layer.

C*********************************************************************
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTAT(MPI_STATUS_SIZE)
      LOGICAL FLAG

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE EXMLST, OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I)

C  LOOP OVER TEMPLATES AND FAULT BLOCKS

      IERR=0
      CALL TIMON(2)
      DO 1 NT=1,NUMTMP
      DO 1 NB=1,NUMBLK
      NM2=MSGS(NT,NB)

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

C  SEND MESSAGE LISTS (NONBLOCKING)

      DO 2 NM=1,NM2
      NE1=MSGSNDF(NM,NT,NB)
      NE2=MSGSNDL(NM,NT,NB)
      MLEN=(NE2-NE1+1)*2

      CALL MPI_ISEND(MSGIND(1,NE1),MLEN,MPI_INTEGER,
     & MSGTRG(NM,NT,NB),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(NM),IERR)

      IF (IERR.GT.0) GO TO 13
    2 CONTINUE

C  RECEIVE MESSAGE LISTS (BLOCKING, ANY SOURCE)

      DO 3 NM=1,NM2
      MLEN=(2000-LSTIND)*2

      CALL MPI_RECV(MSGIND(1,LSTIND+1),MLEN,MPI_INTEGER,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTAT,IERR)
      CALL MPI_GET_COUNT(ISTAT,MPI_INTEGER,NI,IERR)
      MS=ISTAT(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      MSGSRC(NM,NT,NB)=MS
      MSGRCVF(NM,NT,NB)=LSTIND+1
      LSTIND=LSTIND+NI/2
    3 MSGRCVL(NM,NT,NB)=LSTIND

C  WAIT FOR SEND COMPLETION

C      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      DO 7 NM=1,NM2

    9 CALL MPI_TEST(IREQ(NM),FLAG,ISTAT,IERR)
      IF (IERR.GT.0) GO TO 13
      IF (.NOT.FLAG) GO TO 9

    7 CONTINUE

    1 CONTINUE

C  EXITS

      CALL TIMOFF(2)
      IF (LEVELE) WRITE (*,6) MYPRC
    6 FORMAT('MESSAGE LIST EXCHANGE COMPLETE: MYPRC =',I5)

      IF (LEVELE) THEN
         DO 11 KT=1,NUMTMP
         DO 11 NB=1,NUMBLK
         M2=MSGS(KT,NB)
         DO 11 M=1,M2

         LI1=MSGSNDF(M,KT,NB)
         LI2=MSGSNDL(M,KT,NB)

         WRITE(*,16) MYPRC,LI2-LI1+1,MSGTRG(M,KT,NB),NB,KT
   16    FORMAT(' PROC',I4,' SENDING',I5,' COLUMNS TO PROC',I4,
     &      ', BLOCK',I3,', STENCIL',I2)

         IF (BUGKEY(3)) THEN
            LA=LI1
   14       LB=LA+6
            IF (LB.GT.LI2) LB=LI2
            WRITE(*,10) MYPRC,MSGTRG(M,KT,NB),
     &         (MSGIND(1,L),MSGIND(2,L),L=LA,LB)
   10       FORMAT(' SEND FROM',I3,' TO',I3,' JK =',7(2I3,','))
            LA=LB+1
            IF (LB.LT.LI2) GO TO 14
         ENDIF

         LI1=MSGRCVF(M,KT,NB)
         LI2=MSGRCVL(M,KT,NB)

         WRITE(*,17) MYPRC,LI2-LI1+1,MSGTRG(M,KT,NB),NB,KT
   17    FORMAT(' PROC',I4,' RECEIVING',I5,' COLUMNS FROM PROC',I4,
     &      ', BLOCK',I3,', STENCIL',I2)

         IF (BUGKEY(3)) THEN
            LA=LI1
   15       LB=LA+6
            IF (LB.GT.LI2) LB=LI2
            WRITE(*,12) MYPRC,MSGSRC(M,KT,NB),
     &         (MSGIND(1,L),MSGIND(2,L),L=LA,LB)
   12       FORMAT(' RECV AT',I3,' FROM',I3,' JK =',7(2I3,','))
            LA=LB+1
            IF (LB.LT.LI2) GO TO 15
         ENDIF

   11    CONTINUE
      ENDIF

      RETURN
   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN EXMLST FOR PROC',MYPRC
      CALL KILLPRC(IERR)

      END


C*********************************************************************
      SUBROUTINE R4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      REAL*4 AR(IDIM,JDIM,KDIM,LDIM)

      CALL R4UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END


C*********************************************************************
      SUBROUTINE R4UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 REAL*4 array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, REAL*4)

C*********************************************************************
      PARAMETER (NDB=10000000/4)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      REAL*4 AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE R4UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      DO 10 LL=IV1,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_REAL,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) THEN
         NER1=1
         NER2=NFE2-NFE1
         GO TO 13
      ENDIF
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_REAL,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) THEN
         NER1=2
         NER2=LENM
         GO TO 13
      ENDIF
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      NER1=3
      NER2=M2
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN R4UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2,', TAG',MSGTAG(MTM)

      CALL KILLPRC(IERR)

      END

C*********************************************************************
      SUBROUTINE R8UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      REAL*8 AR(IDIM,JDIM,KDIM,LDIM)

      CALL R8UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END



C*********************************************************************
      SUBROUTINE R8UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 REAL*8 array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, REAL*8)

C*********************************************************************

      PARAMETER (NDB=10000000/8)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      REAL*8 AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE R8UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      DO 10 LL=IV1,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_DOUBLE_PRECISION,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) GO TO 13
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_DOUBLE_PRECISION,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN R8UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2
      CALL KILLPRC(IERR)

      END

C*********************************************************************
      SUBROUTINE I4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      INTEGER AR(IDIM,JDIM,KDIM,LDIM)

      CALL I4UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END

C*********************************************************************
      SUBROUTINE I4UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 INTEGER array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, INTEGER)

C  Note: If IV1 = 0 then keyout is being updated and signs are reversed during
C        update.

C*********************************************************************
      PARAMETER (NDB=10000000/4)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      INTEGER AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE I4UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      IV1A=IV1
      IF (IV1A.LT.1) IV1A=1

      DO 10 LL=IV1A,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
      IF (IV1.EQ.0.AND.B(NFE2).GT.0) B(NFE2)=-B(NFE2)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_INTEGER,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) GO TO 13
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_INTEGER,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN I4UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2
      CALL KILLPRC(IERR)
      END

C*********************************************************************
      SUBROUTINE I2UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      INTEGER AR(IDIM,JDIM,KDIM,LDIM)

      CALL I2UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END


C*********************************************************************
      SUBROUTINE I2UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 INTEGER array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, INTEGER)

C  NOTE:  This routine is identical to I4UPDATE since the CRAY has problems
C         with INTEGER*2.  The special treatment of IV1 in I4UPDATE is not
C         retained in I2UPDATE.  The routine is retained in case someone wants
C         to specialize the IPARS for a machine that has no problems with
C         INTEGER*2

C*********************************************************************

      PARAMETER (NDB=10000000/4)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      INTEGER AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE I2UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      DO 10 LL=IV1,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_INTEGER,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) GO TO 13
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_INTEGER,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN I2UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2
      CALL KILLPRC(IERR)
      END

C*********************************************************************
      SUBROUTINE L4UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      LOGICAL AR(IDIM,JDIM,KDIM,LDIM)

      CALL L4UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END

C*********************************************************************
      SUBROUTINE L4UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 LOGICAL array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, LOGICAL)

C*********************************************************************

      PARAMETER (NDB=10000000/4)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      LOGICAL AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE L4UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      DO 10 LL=IV1,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_LOGICAL,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) GO TO 13
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_LOGICAL,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN L4UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2
      CALL KILLPRC(IERR)

      END

C*********************************************************************
      SUBROUTINE L2UPDATE (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR)
C*********************************************************************

CGUS Something stupid to get around EQUIVALENCE AND MODULE incompatibility

C*********************************************************************
      USE scrat1mod

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'

      LOGICAL AR(IDIM,JDIM,KDIM,LDIM)

      CALL L2UPDATE2(IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,A)

      END


C*********************************************************************
      SUBROUTINE L2UPDATE2 (IDIM,JDIM,KDIM,LDIM,IV1,IV2,NBLK,KTMP,AR,B)
C*********************************************************************

C  Exchanges communication layer data for 1 fault block and 1 LOGICAL array.
C  This routine is called only by UPDATE

C  IDIM,JDIM  = First 4 dimensions of AR.  Set LDIM to 1 if there is
C  KDIM,LDIM    no forth dimension. (input,INTEGER)

C  IV1,IV2 = Range to be updated in the fourth dimension. (input,INTEGER)

C  NBLK = Block number (input,INTEGER)

C  KTMP = Update template (input, INTEGER)
C       = 1 ==> adjacent elements only (7 point template)
C       = 2 ==> cube (27 point template) (also cornerless cube)
C       = 3 ==> "higher order" (27 point template)

C  AR = Array to be UPDATED (input, LOGICAL)

C  NOTE: This routine is identical to L4UPDATE since mpi and some machines
C        (CRAY) do not support LOGICAL*2

C*********************************************************************

      PARAMETER (NDB=10000000/4)
      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'
      INCLUDE 'layout.h'
      INCLUDE 'mpif.h'
      INTEGER IREQ(10),ISTATS(MPI_STATUS_SIZE,10)
      INTEGER ISTATR(MPI_STATUS_SIZE)
      LOGICAL AR(IDIM,JDIM,KDIM,LDIM),B(NDB)
!      EQUIVALENCE (A,B)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE L2UPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)
      M2=MSGS(KTMP,NBLK)

      DO 10 LL=IV1,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (M2.LE.0) GO TO 10

C  NONBLOCKING SENDS

      NFE1=1
      DO 2 M=1,M2
      LI1=MSGSNDF(M,KTMP,NBLK)
      LI2=MSGSNDL(M,KTMP,NBLK)
      NFE2=NFE1
      DO 1 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 1 I=1,IDIM
      B(NFE2)=AR(I,J,K,LL)
    1 NFE2=NFE2+1

      CALL MPI_ISEND(B(NFE1),NFE2-NFE1,MPI_LOGICAL,
     & MSGTRG(M,KTMP,NBLK),MSGTAG(MTM),MPI_COMM_WORLD,IREQ(M),IERR)

      IF (IERR.GT.0) GO TO 13
    2 NFE1=NFE2

C  RECEIVE

      LENM=NDB-NFE1+1
      DO 3 M=1,M2

      CALL MPI_RECV(B(NFE1),LENM,MPI_LOGICAL,
     & MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTATR,IERR)
      MS=ISTATR(MPI_SOURCE)

      IF (IERR.GT.0) GO TO 13
      DO 5 MM=1,M2
      IF (MS.EQ.MSGSRC(MM,KTMP,NBLK)) THEN
         LI1=MSGRCVF(MM,KTMP,NBLK)
         LI2=MSGRCVL(MM,KTMP,NBLK)
         GO TO 6
      ENDIF
    5 CONTINUE
      IERR=1
      GO TO 13

    6 NFE2=NFE1
      DO 4 L=LI1,LI2
      J=MSGIND(1,L)-JOFF
      K=MSGIND(2,L)-KOFF
      DO 4 I=1,IDIM
      AR(I,J,K,LL)=B(NFE2)
    4 NFE2=NFE2+1

    3 CONTINUE

C  WAIT FOR SEND COMPLETION

      CALL MPI_WAITALL (M2,IREQ,ISTATS,IERR)

      IF (IERR.GT.0) THEN
         NER1=4
         NER2=M2
         GO TO 13
      ENDIF

C  TERMINATE 4TH INDEX LOOP

   10 CONTINUE

C  EXITS

      RETURN

   13 WRITE (*,*) ' MESSAGE ERROR',IERR,' IN L2UPDATE FOR PROC',MYPRC,
     & ', BLOCK',NBLK,', LOC',NER1,', DATA',NER2
      CALL KILLPRC(IERR)

      END
C*********************************************************************
      SUBROUTINE NOUPDATE (NBLK,IV1,IV2)
C*********************************************************************

C  Steps MSGTAG if no update data is to be transmitted
C  This routine is called only by UPDATE

C*********************************************************************
      INCLUDE 'control.h'

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE NOUPDATE, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      IV1A=IV1
      IF (IV1A.LT.1) IV1A=1

      DO 1 I=IV1A,IV2

      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

    1 CONTINUE

      END
C***********************************************************************
      SUBROUTINE GETMP (I1G,I2G,ITS,IOFF,J1G,J2G,JTS,JOFF,KG,KOFF,
     & IEX,KEYOUT,IDIM,JDIM,KDIM,NBLK,KARY,KIND,ARRY)
C***********************************************************************

C  COPIES DATA TO A PRINT BUFFER FOR GEAOUT().
C  MULTIPROCESSOR ROUTINE NOT A WORK ROUTINE

C***********************************************************************
      INCLUDE 'mpif.h'
      INCLUDE 'control.h'
      INCLUDE 'output.h'
      INCLUDE 'layout.h'
      REAL*4 ARRY(IDIM,JDIM,KDIM,*)
      INTEGER KEYOUT(IDIM,JDIM,KDIM),ISTAT(MPI_STATUS_SIZE)
      INTEGER LMSG(256)

      MTM=MODACT+1

      IF (LEVELE.AND.BUGKEY(1)) WRITE (NFBUG,100)' PROC',MYPRC,
     & ' ENTERING SUBROUTINE GETMP, BLOCK =',NBLK,
     & ', OLD TAG =',MSGTAG(MTM)
 100  FORMAT(A,I,A,I,A,I)

      N=0
      M=N0MAP(NBLK)+KG*NYMAP(NBLK)
!bw for nodal based array
      IF ((KARY.EQ.2).OR.(KARY.EQ.3)) THEN
         M=N0MAPN(NBLK)+KG*(NYMAP(NBLK)+1)
      ENDIF
!bw
      KL=KG-KOFF
      DO 2 I=1,NUMPRC
    2 LMSG(I)=0

C  EXTRACT DATA INTO BUFFER A ON EACH PROCESSOR

      DO 1 ID=I1G,I2G,ITS
      I=ID
      IF (I+ITS.GT.I2G) I=I2G
      IL=I-IOFF
      DO 1 JD=J1G,J2G,JTS
      J=JD
      IF (J+JTS.GT.J2G) J=J2G
      MP=PRCMAP(M+J)
!bw for nodal based array
      IF ((KARY.EQ.2).OR.(KARY.EQ.3)) THEN
         MP=PRCMAPN(M+J)
      ENDIF
!bw
      IF (MP.LT.0) GO TO 1
      LMSG(MP+1)=LMSG(MP+1)+1
      IF (MP.EQ.MYPRC) THEN
         N=LMSG(MYPRC+1)
         JL=J-JOFF
         GO TO (11,12,13,14),KIND
   11    CALL R4GETM(N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
         GO TO 7
   12    CALL R8GETM(N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
         GO TO 7
   13    CALL I2GETM(N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
         GO TO 7
   14    CALL I4GETM(N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
    7    IF (N.EQ.480) GO TO 8
      ENDIF
    1 CONTINUE

C  SEND BUFFER A TO BUFFER B ON PROCESSOR 0

    8 MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

      IF (N.GT.0) THEN

         IF (MYPRC.EQ.0) THEN
            GO TO (21,22,23,23),KIND
   21       DO 24 I=1,N
   24       PBUF4B(I)=PBUF4A(I)
            GO TO 27
   22       DO 25 I=1,N
   25       PBUF8B(I)=PBUF8A(I)
            GO TO 27
   23       DO 26 I=1,N
   26       IPBUF4B(I)=IPBUF4A(I)
            GO TO 27
         ENDIF

         GO TO (31,32,33,33),KIND
   31 CALL MPI_SEND(PBUF4A,N,MPI_REAL,0,MSGTAG(MTM),
     & MPI_COMM_WORLD,IERR)
         GO TO 3
   32 CALL MPI_SEND(PBUF8A,N,MPI_DOUBLE_PRECISION,0,MSGTAG(MTM),
     & MPI_COMM_WORLD,IERR)
         GO TO 3
   33 CALL MPI_SEND(IPBUF4A,N,MPI_INTEGER,0,MSGTAG(MTM),
     & MPI_COMM_WORLD,IERR)

    3    IF (IERR.GT.0) THEN
            WRITE (*,*) ' SEND FAILURE IN GETMP FOR PROC',MYPRC
            CALL KILLPRC(13)
         ENDIF

      ENDIF

   27 IF (MYPRC.EQ.0) THEN

C  CLEAR BUFFER A ON PROC 0

         LP=0
         GO TO (61,62,63,63),KIND

   61    DO 64 ID=I1G,I2G,ITS
         DO 64 JD=J1G,J2G,JTS
         LP=LP+1
   64    PBUF4A(LP)=0.
         GO TO 67

   62    DO 65 ID=I1G,I2G,ITS
         DO 65 JD=J1G,J2G,JTS
         LP=LP+1
   65    PBUF8A(LP)=0.D0
         GO TO 67

   63    DO 66 ID=I1G,I2G,ITS
         DO 66 JD=J1G,J2G,JTS
         LP=LP+1
   66    IPBUF4A(LP)=0

C  RECEIVE IN BUFFER B ON  PROC 0

   67    DO 4 IP=1,NUMPRC
         IF (LMSG(IP).GT.0) THEN
            NS=0
            IF (IP.GT.1) THEN

               GO TO (41,42,43,43),KIND
   41      CALL MPI_RECV(PBUF4B,480,MPI_REAL,
     &     MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTAT,IERR)
               GO TO 5
   42      CALL MPI_RECV(PBUF8B,480,MPI_DOUBLE_PRECISION,
     &     MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTAT,IERR)
               GO TO 5
   43      CALL MPI_RECV(IPBUF4B,480,MPI_INTEGER,
     &     MPI_ANY_SOURCE,MSGTAG(MTM),MPI_COMM_WORLD,ISTAT,IERR)

    5          NS=ISTAT(MPI_SOURCE)
               IF (IERR.GT.0) THEN
                  WRITE (*,*) ' RECEIVE FAILURE IN GETMP FROM PROC',NS
                  CALL KILLPRC(13)
               ENDIF
            ENDIF

C  COPY FROM BUFFER B TO PRINT OFFSET IN BUFFER A (PROC 0)

            LB=0
            LP=0
            DO 6 ID=I1G,I2G,ITS
            I=ID
            IF (I+ITS.GT.I2G) I=I2G
            IL=I-IOFF
            DO 6 JD=J1G,J2G,JTS
            J=JD
            IF (J+JTS.GT.J2G) J=J2G
            LP=LP+1
            MP=PRCMAP(M+J)
!bw for nodal based array
            IF ((KARY.EQ.2).OR.(KARY.EQ.3)) THEN
               MP=PRCMAPN(M+J)
            ENDIF
!bw
            IF (MP.EQ.NS) THEN
               LB=LB+1
               JL=J-JOFF
               GO TO (51,52,53,53),KIND
   51          IF (LB.GT.480) THEN
                  PBUF4A(LP)=-123.
               ELSE
                  PBUF4A(LP)=PBUF4B(LB)
               ENDIF
               GO TO 6
   52          IF (LB.GT.480) THEN
                  PBUF8A(LP)=-123.D0
               ELSE
                  PBUF8A(LP)=PBUF8B(LB)
               ENDIF
               GO TO 6
   53          IF (LB.GT.480) THEN
                  IPBUF4A(LP)=-123
               ELSE
                  IPBUF4A(LP)=IPBUF4B(LB)
               ENDIF
            ENDIF
    6       CONTINUE
         ENDIF
    4    CONTINUE
      ENDIF

      END
C***********************************************************************
      SUBROUTINE R4GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
C***********************************************************************

C  COPIES A REAL*4 VARIABLE TO A BUFFER FOR GETMP().

C***********************************************************************
      INCLUDE 'output.h'
      REAL*4 ARRY(IDIM,JDIM,KDIM,*)
      INTEGER KEYOUT(IDIM,JDIM,KDIM)

      IF (KARY.EQ.2.OR.KEYOUT(IL,JL,KL).EQ.1) THEN
         PBUF4A(N)=ARRY(IL,JL,KL,IEX)
      ELSE
         PBUF4A(N)=0.
      ENDIF

      END
C***********************************************************************
      SUBROUTINE R8GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
C***********************************************************************

C  COPIES A REAL*8 VARIABLE TO A BUFFER FOR GETMP().

C***********************************************************************
      INCLUDE 'output.h'
      REAL*8 ARRY(IDIM,JDIM,KDIM,*)
      INTEGER KEYOUT(IDIM,JDIM,KDIM)

      IF (KARY.EQ.2.OR.KEYOUT(IL,JL,KL).EQ.1) THEN
         PBUF8A(N)=ARRY(IL,JL,KL,IEX)
      ELSE
         PBUF8A(N)=0.D0
      ENDIF

      END
C***********************************************************************
      SUBROUTINE I2GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
C***********************************************************************

C  COPIES A INTEGER VARIABLE TO A BUFFER FOR GETMP().

C***********************************************************************
      INCLUDE 'output.h'
      INTEGER ARRY(IDIM,JDIM,KDIM,*)
      INTEGER KEYOUT(IDIM,JDIM,KDIM)

      IF (KARY.EQ.2.OR.KEYOUT(IL,JL,KL).EQ.1) THEN
         IPBUF4A(N)=ARRY(IL,JL,KL,IEX)
      ELSE
         IPBUF4A(N)=0
      ENDIF

      END
C***********************************************************************
      SUBROUTINE I4GETM (N,IL,JL,KL,IDIM,JDIM,KDIM,KEYOUT,ARRY,KARY,IEX)
C***********************************************************************

C  COPIES A INTEGER VARIABLE TO A BUFFER FOR GETMP().

C***********************************************************************
      INCLUDE 'output.h'
      INTEGER ARRY(IDIM,JDIM,KDIM,*)
      INTEGER KEYOUT(IDIM,JDIM,KDIM)

      IF (KARY.EQ.2.OR.KEYOUT(IL,JL,KL).EQ.1) THEN
         IPBUF4A(N)=ARRY(IL,JL,KL,IEX)
      ELSE
         IPBUF4A(N)=0
      ENDIF

      END
