C  INFCOMM.F - DUAL APPROXIMATION MULTIBLOCK COMMUNICATIONS AND UTILITIES

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE POLLY    (AX,BX,BX,BY,NP,PX,PY,TOL)
C  SUBROUTINE SORTIF   ()
C  SUBROUTINE IFCOMM   (NERR)
C  SUBROUTINE IFEXCD   (NERR)
C  SUBROUTINE RIDAT    (IB,NERR)
C  SUBROUTINE IFGES8   (N_A,N_BUF,NERR)
C  SUBROUTINE IFGES4   (N_A,N_BUF,NERR)
C  SUBROUTINE IFTOBUF8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                       KL2,KEYOUT,NBLK,A,B)
C  SUBROUTINE IFTOBUF4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                       KL2,KEYOUT,NBLK,A,B)
C  SUBROUTINE IFTRAN   ()

C  CODE HISTORY:

C  JOHN WHEELER     2/19/99     ALPHA CODE
C  JOHN WHEELER     3/14/99     MEMORY MANAGEMENT FOR INTERFACE BUFFERS
C  SUNIL G THOMAS   9/--/09     ADDED UTILITY TO ENABLE VELOCITY VIS
C                               WHEN EVMFEM IS ACTIVE
C  BEN GANIS         12/2/15    LOCAL FLUX IMPLEMENTATION FOR MULTIBLOCK
C                               WITH GENERAL HEXAHEDRA

C*********************************************************************
      SUBROUTINE POLLY (AX,AY,BX,BY,NP,PX,PY,TOL)
C*********************************************************************

C  DETERMINES THE POLYGON OF OVERLAP FOR TWO RECTANGLES

C  AX(),AY() = CORNERS OF A RECTANGLE WITH SIDES PARALLEL TO THE COORDINATE
C              AXISES.  THE POINTS ARE ORDERED IN CLOCKWISE SEQUENCE WITH
C              THE FIRST SIDE PARALLEL TO THE Y AXIS. (INPUT, REAL*8)

C  BX(),BY() = CORNERS OF A RECTANGLE WITH POINTS IN CLOCKWISE SEQUENCE.
C              (INPUT, REAL*8).  BX AND BY MAY BE MODIFIED TO ELIMINATE
C              SMALL DIFFERENCES BETWEEN RECTANGLES A AND B.

C  NP = NUMBER OF CORNERS IN THE OVERLAP POLYGON. (OUTPUT, INTEGER)

C  PX(),PY() = CORNERS IN THE OVERLAP POLYGON. (OUTPUT, REAL*8)
C              PX(NP+1),PY(NP+1) WILL EQUAL PX(1),PY(2)

C  TOL = TOLLERANCE BELOW WHICH COORDINATE DIFFERENCES MAY BE NEGLECTED

C*********************************************************************
      IMPLICIT NONE
      REAL*8 AX(4),AY(4),BX(4),BY(4),PX(9),PY(9),TOL,X,Y,D,DD
      REAL*8 T1,T2,T3,T4
      INTEGER IA(4),IB(4),I,K,JA,JB,NA,NB,NC,NP,NAP,NBP

C  ELIMINATE NEAR MISSES

      DO 1 JB=1,4
      IF (ABS(BX(JB)-AX(1)).LE.TOL) BX(JB)=AX(1)
      IF (ABS(BX(JB)-AX(3)).LE.TOL) BX(JB)=AX(3)
      IF (ABS(BY(JB)-AY(1)).LE.TOL) BY(JB)=AY(1)
      IF (ABS(BY(JB)-AY(3)).LE.TOL) BY(JB)=AY(3)
    1 CONTINUE

C  LOOK FOR B VERTICES CONTAINED IN A

      DO 2 JB=1,4
      IB(JB)=0
      IF ((BX(JB)-AX(1))*(BX(JB)-AX(3)).LE.0.D0.AND.
     & (BY(JB)-AY(2))*(BY(JB)-AY(4)).LE.0.D0) IB(JB)=1
    2 CONTINUE

C  LOOK FOR A VERTICES CONTAINED IN B

      DO 3 JA=1,4
      IA(JA)=0
      X=AX(JA)
      Y=AY(JA)
      T1=(X-BX(1))*(BY(2)-BY(1))-(Y-BY(1))*(BX(2)-BX(1))
      T3=(X-BX(3))*(BY(4)-BY(3))-(Y-BY(3))*(BX(4)-BX(3))
      IF (T1*T3.GE.0.D0) THEN
         T2=(X-BX(2))*(BY(3)-BY(2))-(Y-BY(2))*(BX(3)-BX(2))
         T4=(X-BX(4))*(BY(1)-BY(4))-(Y-BY(4))*(BX(1)-BX(4))
         IF (T2*T4.GE.0.D0) THEN
            DO 4 JB=1,4
            IF (X.EQ.BX(JB).AND.Y.EQ.BY(JB)) GO TO 3
    4       CONTINUE
            IA(JA)=1
         ENDIF
      ENDIF
    3 CONTINUE

C  LOOP TO FIND OVERLAP

      NP=0
      NA=1
      NB=1
      K=0
      NC=0

    5 NAP=MOD(NA,4)+1
      NBP=MOD(NB,4)+1

C  LOOK FOR INTERSECTIONS

      I=0
      IF (NA.EQ.1.OR.NA.EQ.3) THEN
         X=AX(NA)
         D=X-BX(NB)
         IF (D*(X-BX(NBP)).LT.0.D0) THEN
            DD=BX(NBP)-BX(NB)
            IF (ABS(DD).GT.TOL) THEN
               Y=BY(NB)+D*(BY(NBP)-BY(NB))/DD
               IF ((Y-AY(NA))*(Y-AY(NAP)).LT.0) I=1
            ENDIF
         ENDIF
      ELSE
         Y=AY(NA)
         D=Y-BY(NB)
         IF (D*(Y-BY(NBP)).LT.0.D0) THEN
            DD=BY(NBP)-BY(NB)
            IF (ABS(DD).GT.TOL) THEN
               X=BX(NB)+D*(BX(NBP)-BX(NB))/DD
               IF ((X-AX(NA))*(X-AX(NAP)).LT.0) I=1
            ENDIF
         ENDIF
      ENDIF
      IF (NA.EQ.2.OR.NA.EQ.3) D=-D

C  RECORD AND ADVANCE

      IF (I.EQ.1) THEN
         NP=NP+1
         PX(NP)=X
         PY(NP)=Y
         IF (NP.GT.1) THEN
           IF (PX(NP).EQ.PX(1).AND.PY(NP).EQ.PY(1)) GO TO 6
         ENDIF
         IF (D.GT.0.D0) THEN
            NA=NAP
            K=-1
         ELSE
            NB=NBP
            K=1
         ENDIF
      ENDIF
      IF (K.GE.0.AND.IA(NAP).EQ.1) THEN
         NP=NP+1
         PX(NP)=AX(NAP)
         PY(NP)=AY(NAP)
         NA=NAP
         K=1
      ENDIF
      IF (K.LE.0.AND.IB(NBP).EQ.1) THEN
         NP=NP+1
         PX(NP)=BX(NBP)
         PY(NP)=BY(NBP)
         NB=NBP
         K=-1
      ENDIF

      IF (NP.GT.1) THEN
        IF (PX(NP).EQ.PX(1).AND.PY(NP).EQ.PY(1)) GO TO 6
      ENDIF
      NC=NC+1
      IF (NC.GT.15) GO TO 6
      IF (NA.NE.NAP.AND.NB.NE.NBP) THEN
         IF (K.EQ.0) THEN
            IF (NBP.EQ.1) THEN
               IF (NAP.EQ.4) GO TO 6
               NA=NAP
            ENDIF
            NB=NBP
         ELSE
            IF (K.GT.0) THEN
               IF (IB(NBP).EQ.1) THEN
                  K=-1
               ELSE
                  NB=NBP
               ENDIF
            ELSE
               IF (IA(NAP).EQ.1) THEN
                  K=1
               ELSE
                  NA=NAP
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      GO TO 5

C  RESULT

    6 NP=NP-1
      IF (NP.LT.3) NP=0

      END
C*********************************************************************
      SUBROUTINE SORTIF()
C*********************************************************************

C  Sorts Interface data by source element and completes definition of
C  interface base data

C*********************************************************************
      USE dualmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'

      INCLUDE 'sblkc.h'

      INTEGER I,J,K,K2,N1,N2,I1,NB
      INTEGER TMP(2,4,3),TMP2(4)      ! bag8
      REAL*8 DUM

C  SORT TFINS, KDIRS, IJKT, JBLOCK, AND NAEB, BY VALUES IN NAEB TO GROUP
C  DATA RELATED TO THE SAME A BLOCK ELEMENT (SAME k VALUE)

      DO 5 NB=1,NUMBLK
      IF (NIEBS(NB).LT.2) GO TO 5

      K=IIEBS(NB)
      K2=K+NIEBS(NB)-1
      N1=ICGES(K)
      N2=NFICGES-1

    1 I1=N1+1
      DO 2 I=I1,N2
      IF (NAEB(I).EQ.K) THEN
         N1=N1+1
         IF (I.GT.N1) THEN

            NAEB(I)=NAEB(N1)
            NAEB(N1)=K

            DUM=TFINS(I)
            TFINS(I)=TFINS(N1)
            TFINS(N1)=DUM

C COPY FOR AREAS
            DUM=AREAI(I)
            AREAI(I)=AREAI(N1)
            AREAI(N1)=DUM

C COPY FOR ISOLATED MATRIX SOLVES (NON-MODEL DEPENDENT)
            DUM=TDDFINS0(I)
            TDDFINS0(I)=TDDFINS0(N1)
            TDDFINS0(N1)=DUM
C

            J=IJKT(1,I)
            IJKT(1,I)=IJKT(1,N1)
            IJKT(1,N1)=J
            J=IJKT(2,I)
            IJKT(2,I)=IJKT(2,N1)
            IJKT(2,N1)=J
            J=IJKT(3,I)
            IJKT(3,I)=IJKT(3,N1)
            IJKT(3,N1)=J

C COPY FOR ISOLATED MATRIX SOLVES (NON-MODEL DEPENDENT)
            J=TRDDIJKT(1,I)
            TRDDIJKT(1,I)=TRDDIJKT(1,N1)
            TRDDIJKT(1,N1)=J
            J=TRDDIJKT(2,I)
            TRDDIJKT(2,I)=TRDDIJKT(2,N1)
            TRDDIJKT(2,N1)=J
            J=TRDDIJKT(3,I)
            TRDDIJKT(3,I)=TRDDIJKT(3,N1)
            TRDDIJKT(3,N1)=J
C

            J=KDIRS(I)
            KDIRS(I)=KDIRS(N1)
            KDIRS(N1)=J

            J=JBLOCK(I)
            JBLOCK(I)=JBLOCK(N1)
            JBLOCK(N1)=J

! bag8 - also sort DFAC2IJK and DFAC2NOD
            IF (ALLOCATED(DFAC2IJK)) THEN
              TMP(:,:,:)=DFAC2IJK(I,:,:,:)
              DFAC2IJK(I,:,:,:)=DFAC2IJK(N1,:,:,:)
              DFAC2IJK(N1,:,:,:)=TMP(:,:,:)
              TMP2(:)=DFAC2NOD(I,:)
              DFAC2NOD(I,:)=DFAC2NOD(N1,:)
              DFAC2NOD(N1,:)=TMP2(:)
            ENDIF

         ENDIF
      ENDIF

    2 CONTINUE
      IF (N1.LT.N2-2) THEN
         N1=N1+1
         K=NAEB(N1)
         GO TO 1
      ENDIF

C  RESET ICGES

      K=IIEBS(NB)
      N1=ICGES(K)+1
      DO 4 I=N1,N2
      IF (NAEB(I).NE.NAEB(I-1)) ICGES(NAEB(I))=I
    4 CONTINUE

    5 CONTINUE

      END
C*********************************************************************
      SUBROUTINE IFCOMM (NERR)
C*********************************************************************

C  Sets up interface communications
C  Note that a processor must send messages to itself if it has grid
C  elements on both sides of an interface

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER J,K,M,JB,KB,NA,NB,NPB,J1,J2,K1,K2,MP,NERR,L,L1,L2,NBJ,NPJ
      NFESR=1

C  COUNT THE NUMBER OF MESSAGES AND IDENTIFY TARGET PROCESSORS

      DO 1 NA=1,NUMBLK
      NPAI(NA)=0
      IF (NIEBS(NA).EQ.0) GO TO 1
      K1=IIEBS(NA)
      K2=K1+NIEBS(NA)-1
      DO 2 K=K1,K2
      J1=ICGES(K)
      J2=J1+NCGES(K)-1
      DO 3 J=J1,J2
      NB=JBLOCK(J)
      JB=IJKT(2,J)
      KB=IJKT(3,J)
      NPB=PRCMAP(N0MAP(NB)+KB*NYMAP(NB)+JB)
      DO 4 M=1,NPAI(NA)
      IF (NPB.EQ.NPSRI(M,NA).AND.NB.EQ.NBSRI(M,NA)) GO TO 3
    4 CONTINUE
      MP=NPAI(NA)+1
      IF (MP.GT.10) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,10)
   10    FORMAT (/' ERROR # 508; MAX NUMBER OF INTERFACE MESSAGES',
     &   ' PER BLOCK EXCEEDED')
         RETURN
      ENDIF
      NPAI(NA)=MP
      NBSRI(MP,NA)=NB
      NPSRI(MP,NA)=NPB
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE

C  IDENTIFY ELEMENTS IN EACH MESSAGE

      DO 8 NA=1,NUMBLK
      IF (NIEBS(NA).EQ.0) GO TO 8
      K1=IIEBS(NA)
      K2=K1+NIEBS(NA)-1
      DO 5 M=1,NPAI(NA)
      NB=NBSRI(M,NA)
      NPB=NPSRI(M,NA)
      NESNDI(M,NA)=0
      IESNDI(M,NA)=NFESR
      L1=NFESR
      DO 5 K=K1,K2
      J1=ICGES(K)
      J2=J1+NCGES(K)-1
      DO 6 J=J1,J2
      NBJ=JBLOCK(J)
      JB=IJKT(2,J)
      KB=IJKT(3,J)
      NPJ=PRCMAP(N0MAP(NBJ)+KB*NYMAP(NBJ)+JB)
      IF (NPJ.NE.NPB.OR.NBJ.NE.NB) GO TO 6
      L2=L1+NESNDI(M,NA)-1
      DO 7 L=L1,L2
      IF (K.EQ.KFESR(L)) GO TO 6
    7 CONTINUE
      IF (NFESR.GT.90000) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,9)
    9    FORMAT (/' ERROR # 507; MAX NUMBER OF INTERFACE MESSAGE',
     &   ' ENTRIES EXCEEDED')
         RETURN
      ENDIF
      KFESR(NFESR)=K
      NFESR=NFESR+1
      NESNDI(M,NA)=NESNDI(M,NA)+1
    6 CONTINUE
    5 CONTINUE
    8 CONTINUE

C  INFORM RECEVING PROCESSORS OF ELEMENT SEQUENCES IN MESSAGES

      CALL IFEXCD (NERR)

      END
C*********************************************************************
      SUBROUTINE IFEXCD (NERR)
C*********************************************************************

C  Exchanges interface data between processors

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER J,K,L,L1,L2,M,MTM,MM,NA,IOFFA,JOFFA,KOFFA,MERR,NERR
      INTEGER IB(3,90000+2)

      MTM=MODACT+1
      MSGTAG(MTM)=MSGTAG(MTM)+1
      IF (MSGTAG(MTM).GT.MSGTAG2(MTM)) MSGTAG(MTM)=MSGTAG1(MTM)

C  LOOP OVER THE FAULT BLOCKS AND SEND MESSAGES

      DO 1 NA=1,NUMBLK
      IF (NIEBS(NA).EQ.0) GO TO 1
      CALL BLKOFF(NA,IOFFA,JOFFA,KOFFA,MERR)
      DO 2 M=1,NPAI(NA)

C  LOAD SOURCE BLOCK NUMBER, TARGET BLOCK NUMBER, NUMBER OF MESSAGE ENTRIES,
C  SOURCE PROCESSOR, TARGET PROCESSOR, SOURCE OFFSET AND IJKs INTO IB

      IB(1,1)=NA
      IB(2,1)=NBSRI(M,NA)
      IB(3,1)=NESNDI(M,NA)
      IB(1,2)=MYPRC
      IB(2,2)=NPSRI(M,NA)
      IB(3,2)=IESNDI(M,NA)

      L1=IESNDI(M,NA)
      L2=L1+NESNDI(M,NA)-1
      J=2
      DO 3 L=L1,L2
      J=J+1
      K=KFESR(L)
      IB(1,J)=IJKS(1,K)+IOFFA
      IB(2,J)=IJKS(2,K)+JOFFA
    3 IB(3,J)=IJKS(3,K)+KOFFA

C  SEND TO ANOTHER PROCESSOR OR PROCESS ON THE CURRENT PROCESSOR

      IF (MYPRC.EQ.NPSRI(M,NA)) THEN
         CALL RIDAT (IB,NERR)
      ELSE
      CALL XPIDAT (IB,NERR)
      ENDIF

    2 CONTINUE
    1 CONTINUE

C  RECEIVE AND PROCESS DATA FROM OTHER PROCESSORS (IN ANY ORDER)

      DO 4 NA=1,NUMBLK
      IF (NIEBS(NA).EQ.0) GO TO 4
      DO 5 MM=1,NPAI(NA)
      IF (MYPRC.EQ.NPSRI(MM,NA)) GO TO 5
      CALL IPIDAT (IB,NERR)
      CALL RIDAT (IB,NERR)
      IF (NERR.GT.0) RETURN
    5 CONTINUE
    4 CONTINUE

      END
C*********************************************************************
      SUBROUTINE RIDAT (IB,NERR)
C*********************************************************************

C  PROCESSES INTERFACE INITIALIZATION DATA (ASSIGN LIBUF())

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C  IB() = DATA BUFFER (INPUT, INTEGER)

C  IB(1,1)=SOURCE BLOCK
C  IB(2,1)=TARGET BLOCK
C  IB(3,1)=NUMBER ELEMENTS IN THE MESSAGE
C  IB(1,2)=SOURCE PROCESSOR
C  IB(2,2)=TARGET PROCESSOR
C  IB(3,2)=SOURCE OFFSET OF FIRST ELEMENT
C  IB( ,n+2) = IJK OF MESSAGE ELEMENT n IN TARGET BLOCK (GLOBAL)

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER NBR,NBS,NPR,NPS,NSO,I,II,J,J1,J2,K,K1,K2,L,LI,L2,NERR
      INTEGER IB(3,90000+2)

      NBS=IB(1,1)
      NBR=IB(2,1)
      L2=IB(3,1)+2
      NPS=IB(1,2)
      NPR=IB(2,2)
      NSO=IB(3,2)

C  IDENTIFY MESSAGE NUMBER AND INITIALIZE OFFSET

      DO 5 II=1,NPAI(NBR)
      I=II
C      IF (NPS.EQ.NPSRI(II,NBR)) GO TO 6
      IF (NPS.EQ.NPSRI(II,NBR) .AND. NBS.EQ.NBSRI(II,NBR)) GO TO 6
    5 CONTINUE

      NERR=NERR+1
      WRITE (*,7) MYPRC,NBS,NBR,NPS,NPR
      IF (LEVELC) WRITE (NFOUT,7) MYPRC,NBS,NBR,NPS,NPR
    7 FORMAT (/' ERROR # 509; UNEXPECTED INTERFACE MESSAGE',5I5)
      RETURN

    6 IF (NPS.EQ.MYPRC) THEN
         LI=NSO
      ELSE
         LI=NFESR
         NFESR=NFESR+IB(3,1)
      ENDIF
      NERECI(I,NBR)=IB(3,1)
      IERECI(I,NBR)=LI

C  SET LIBUF()

      K1=IIEBS(NBR)
      K2=K1+NIEBS(NBR)-1
      DO 1 K=K1,K2
      J1=ICGES(K)
      J2=J1+NCGES(K)-1
      DO 2 J=J1,J2
      IF (JBLOCK(J).NE.NBS) GO TO 2
      DO 3 L=3,L2
      IF (IB(1,L).EQ.IJKT(1,J).AND.IB(2,L).EQ.IJKT(2,J).AND.
     &   IB(3,L).EQ.IJKT(3,J)) LIBUF(J)=LI+L-3
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE

      END
C*********************************************************************
      SUBROUTINE IFGES8 (N_A,NBUF,NERR)
C*********************************************************************

C  Loads an interface buffer from a REAL*8 grid-element array.  Does not
C  (repeat not) send it to other processors.  This routine can (and should)
C  be used even on single processor machines.  The routine can be called
C  only from executive routines.  Interfaces on all fault blocks are
C  processed by a single call.

C  N_A = Array number of the REAL*8 grid element array that is the source
C        of the transfer (input, INTEGER).

C  NBUF = Buffer number of a REAL*8 buffer used in the transfer
C         (input, INTEGER).

C  NERR = Error number stepped by 1 on error (input & output, INTEGER)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'sblkc.h'

      INTEGER N_A,NBUF,NERR,IA(4)

      EXTERNAL IFTOBUF8

C  LOAD BUFFER VIA A WORK ROUTINE

      IA(1)=3
      IA(2)=N_A
      IA(3)=N_BUFDIM
      IA(4)=N_BUFIF
      NUMBUFU=NBUF

      CALL CALLWORK (IFTOBUF8,IA)

      END

C*********************************************************************
      SUBROUTINE IFGES4 (N_A,NBUF,NERR)
C*********************************************************************

C  Loads an interface buffer from a REAL*4 grid-element array.  Does not
C  (repeat not) send it to other processors.  This routine can (and should)
C  be used even on single processor machines.  The routine can be called
C  only from executive routines.  Interfaces on all fault blocks are
C  processed by a single call.

C  N_A = Array number of the REAL*4 grid ellement array that is the source
C        of the transfer (input, INTEGER).

C  NBUF = Buffer number of a REAL*4 buffer used in the transfer
C         (input, INTEGER).

C  NERR = Error number stepped by 1 on error (input & output, INTEGER)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'sblkc.h'

      INTEGER N_A,NBUF,NERR,IA(4)

      EXTERNAL IFTOBUF4

C  LOAD BUFFER VIA A WORK ROUTINE

      IA(1)=3
      IA(2)=N_A
      IA(3)=N_BUFDIM
      IA(4)=N_BUFIF
      NUMBUFU=NBUF

      CALL CALLWORK (IFTOBUF4,IA)

      END
C*********************************************************************
      SUBROUTINE IFTOBUF8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &           KL2,KEYOUT,NBLK,A,NBUFDIM,BUFIF8)
C*********************************************************************

C  Copies interface data from a REAL*8 grid-element array to a buffer.
C  This is a work routine.

C  A() = Grid-element array (input, REAL*8)

C  NBUFDIM = First dimension of BUFIF8() (input, INTEGER)

C  BUFIF8(,) = Buffer (output, REAL*8)

C  NOTE: The buffer number is passed in NUMBUFU

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'sblkc.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,NBUFDIM
      INTEGER JL1V(KDIM), JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      INTEGER II,II1,II2,K,M
      REAL*8 A(IDIM,JDIM,KDIM), BUFIF8(NBUFDIM,*)

      IF (NIEBS(NBLK).EQ.0) RETURN

      DO 1 M=1,NPAI(NBLK)
      II1=IESNDI(M,NBLK)
      II2=II1+NESNDI(M,NBLK)-1
      DO 1 II=II1,II2
      K=KFESR(II)
    1 BUFIF8(II,NUMBUFU)=A(IJKS(1,K),IJKS(2,K),IJKS(3,K))

      END
C*********************************************************************
      SUBROUTINE IFTOBUF4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &           KL2,KEYOUT,NBLK,A,NBUFDIM,BUFIF4)
C*********************************************************************

C  Copies interface data from a REAL*8 grid-element array to a buffer.
C  This is a work routine.

C  A() = Grid-element array (input, REAL*4)

C  NBUFDIM = 1/2 the first dimension of BUFIF4() (input, INTEGER)

C  BUFIF4(,) = Buffer (output, REAL*4)

C  NOTE: The buffer number is passed in NUMBUFU

C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'sblkc.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,NBUFDIM
      INTEGER JL1V(KDIM), JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      INTEGER M,NES,II,II1,II2,K
CSGT      DIMENSION A(IDIM,JDIM,KDIM), BUFIF4(2*NBUFDIM,*)
      REAL*4 A(IDIM,JDIM,KDIM), BUFIF4(2*NBUFDIM,*)

      IF (NIEBS(NBLK).EQ.0) RETURN

      DO 1 M=1,NPAI(NBLK)
      NES=NESNDI(M,NBLK)
      II1=IESNDI(M,NBLK)
      II2=II1+NES-1
      DO 1 II=II1,II2
      K=KFESR(II)
    1 BUFIF4(II,NUMBUFU)=A(IJKS(1,K),IJKS(2,K),IJKS(3,K))

      END

C*********************************************************************
      SUBROUTINE IFTRAN ()
C*********************************************************************

C  Completes definition of interface transmissability.
C  Computes depth differences between interface elements

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER I,IA(6),NERR,NBEM(19)
      EXTERNAL IFTRANW

C  PUT PERMEABILITIES FOR THE B BLOCK(S) IN BUFFERS 1 TO 3

      NERR=0
      CALL IFGES4 (N_XPERM,1,NERR)
      CALL IFGES4 (N_YPERM,2,NERR)
      CALL IFGES4 (N_ZPERM,3,NERR)

      DO 1 I=1,19
    1 NBEM(I)=3

      CALL PIFBUF4(NBEM,NERR)

      IF (NERR.GT.0) RETURN

C  COMPUTE TRANSMISSABILITIES AND DEPTH CHANGES IN A WORK ROUTINE

      IA(1)=5
      IA(2)=N_XPERM
      IA(3)=N_YPERM
      IA(4)=N_ZPERM
      IA(5)=N_BUFDIM
      IA(6)=N_BUFIF
      CALL CALLWORK(IFTRANW,IA)

      END
C*********************************************************************
      SUBROUTINE IFTRANW (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &           KL2,KEYOUT,NBLK,PX,PY,PZ,NBUFDIM,BUFIF4)
C*********************************************************************

C  Work routine for computing interface transmissabilities.

C  PX() = Permeabilities (input, REAL*4)
C  PY()
C  PZ()

C  TX() = Transmissabilities (output, REAL*8)
C  TY()
C  TZ()

C  NBUFDIM = 1/2 the first dimension of BUFIF4(,) (input, INTEGER)

C  BUFIF4(i,j) = Location i in interface buffer j (input, REAL*4)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,NBUFDIM
      INTEGER JL1V(KDIM), JL2V(KDIM), KEYOUT(IDIM,JDIM,KDIM)
      INTEGER IOFFA,JOFFA,KOFFA,IA,JA,KA,IAG,JAG,KAG,IBG,JBG,KBG,
     &  K,K1,K2,J,J1,J2,MERR,IB,NB
      REAL*8 PA,PB,DA,DB,CVC
      REAL*4 PX(IDIM,JDIM,KDIM),PY(IDIM,JDIM,KDIM),PZ(IDIM,JDIM,KDIM),
     &       BUFIF4(2*NBUFDIM,*)

c bag8
      REAL*8 XA(3,8),XB(3,8),XTMP,YTMP,ZTMP,NMLA(3),NMLB(3)
      INTEGER KDB

C  CVC = 2 * 24 * 2.636786E-4 ==> TRAN = cu-ft cp / psi day

C bag8
C      DATA CVC/.12656573D-1/
      CVC = 2 * CONV_FACTOR

      IF (NIEBS(NBLK).EQ.0) RETURN

      CALL BLKOFF(NBLK,IOFFA,JOFFA,KOFFA,MERR)

      K1=IIEBS(NBLK)
      K2=K1+NIEBS(NBLK)-1
      DO 1 K=K1,K2
      IA=IJKS(1,K)
      JA=IJKS(2,K)
      KA=IJKS(3,K)
      IAG=IA+IOFFA
      JAG=JA+JOFFA
      KAG=KA+KOFFA

      J1=ICGES(K)
      J2=J1+NCGES(K)-1
      DO 1 J=J1,J2
      NB=JBLOCK(J)
      IBG=IJKT(1,J)
      JBG=IJKT(2,J)
      KBG=IJKT(3,J)
      IB=LIBUF(J)

C bag8
C      IF (EVFEM_HEX.GE.2) THEN

C       XTMP=XREC(IAG,NBLK);YTMP=YREC(JAG,NBLK);ZTMP=ZREC(KAG,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,1),XA(2,1),XA(3,1))
C       XTMP=XREC(IAG+1,NBLK);YTMP=YREC(JAG,NBLK);ZTMP=ZREC(KAG,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,2),XA(2,2),XA(3,2))
C       XTMP=XREC(IAG+1,NBLK);YTMP=YREC(JAG+1,NBLK);ZTMP=ZREC(KAG,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,3),XA(2,3),XA(3,3))
C       XTMP=XREC(IAG,NBLK);YTMP=YREC(JAG+1,NBLK);ZTMP=ZREC(KAG,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,4),XA(2,4),XA(3,4))
C       XTMP=XREC(IAG,NBLK);YTMP=YREC(JAG,NBLK);ZTMP=ZREC(KAG+1,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,5),XA(2,5),XA(3,5))
C       XTMP=XREC(IAG+1,NBLK);YTMP=YREC(JAG,NBLK);ZTMP=ZREC(KAG+1,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,6),XA(2,6),XA(3,6))
C       XTMP=XREC(IAG+1,NBLK);YTMP=YREC(JAG+1,NBLK);ZTMP=ZREC(KAG+1,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,7),XA(2,7),XA(3,7))
C       XTMP=XREC(IAG,NBLK);YTMP=YREC(JAG+1,NBLK);ZTMP=ZREC(KAG+1,NBLK)
C       CALL MAP(XTMP,YTMP,ZTMP,XA(1,8),XA(2,8),XA(3,8))
C       CALL COMP_NORMAL_AND_DX(XA,KDIRS(J),NMLA,DA)
C       PA=NMLA(1)*PX(IA,JA,KA)+NMLA(2)*PY(IA,JA,KA)+
C     &    NMLA(3)*PZ(IA,JA,KA)

C       XTMP=XREC(IBG,NB);YTMP=YREC(JBG,NB);ZTMP=ZREC(KBG,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,1),XB(2,1),XB(3,1))
C       XTMP=XREC(IBG+1,NB);YTMP=YREC(JBG,NB);ZTMP=ZREC(KBG,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,2),XB(2,2),XB(3,2))
C       XTMP=XREC(IBG+1,NB);YTMP=YREC(JBG+1,NB);ZTMP=ZREC(KBG,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,3),XB(2,3),XB(3,3))
C       XTMP=XREC(IBG,NB);YTMP=YREC(JBG+1,NB);ZTMP=ZREC(KBG,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,4),XB(2,4),XB(3,4))
C       XTMP=XREC(IBG,NB);YTMP=YREC(JBG,NB);ZTMP=ZREC(KBG+1,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,5),XB(2,5),XB(3,5))
C       XTMP=XREC(IBG+1,NB);YTMP=YREC(JBG,NB);ZTMP=ZREC(KBG+1,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,6),XB(2,6),XB(3,6))
C       XTMP=XREC(IBG+1,NB);YTMP=YREC(JBG+1,NB);ZTMP=ZREC(KBG+1,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,7),XB(2,7),XB(3,7))
C       XTMP=XREC(IBG,NB);YTMP=YREC(JBG+1,NB);ZTMP=ZREC(KBG+1,NB)
C       CALL MAP(XTMP,YTMP,ZTMP,XB(1,8),XB(2,8),XB(3,8))
C       IF (KDIRS(J).EQ.1) THEN; KDB=4
C       ELSEIF (KDIRS(J).EQ.2) THEN; KDB=5
C       ELSEIF (KDIRS(J).EQ.3) THEN; KDB=6
C       ELSEIF (KDIRS(J).EQ.4) THEN; KDB=1
C       ELSEIF (KDIRS(J).EQ.5) THEN; KDB=2
C       ELSEIF (KDIRS(J).EQ.6) THEN; KDB=3
C       ENDIF
C       CALL COMP_NORMAL_AND_DX(XB,KDB,NMLB,DB)
C       PB=NMLB(1)*BUFIF4(IB,1)+NMLB(2)*BUFIF4(IB,2)+
C     &    NMLB(3)*BUFIF4(IB,3)

C       GOTO 3

C      ELSE

        GO TO (11,12,13,11,12,13),KDIRS(J)

C      ENDIF

   11 PA=PX(IA,JA,KA)
      DA=DXREC(IAG,NBLK)
      PB=BUFIF4(IB,1)
      DB=DXREC(IBG,NB)
      GO TO 3

   12 PA=PY(IA,JA,KA)
      DA=DYREC(JAG,NBLK)
      PB=BUFIF4(IB,2)
      DB=DYREC(JBG,NB)
      GO TO 3

   13 PA=PZ(IA,JA,KA)
      DA=DZREC(KAG,NBLK)
      PB=BUFIF4(IB,3)
      DB=DZREC(KBG,NB)

    3 IF (PA.GT.0.D0.AND.PB.GT.0.D0) THEN
!         TFINS(J)=CVC*TFINS(J)/(DA/PA+DB/PB)
         TFINS(J)=CVC*AREAI(J)/(DA/PA+DB/PB)    ! bag8
      ELSE
         TFINS(J)=0.D0
      ENDIF

    1 CONTINUE

      END

!----------------------------------------------------------------------
! bag8
!   Hard coded mapping for first EvFEM test cases.  This will be
!   replaced later.
!
      SUBROUTINE MAP(XHAT,YHAT,ZHAT,X,Y,Z)
      IMPLICIT NONE
      INCLUDE 'sblkc.h'
      REAL*8, INTENT(IN) :: XHAT,YHAT,ZHAT
      REAL*8, INTENT(OUT) :: X,Y,Z
      REAL*8, PARAMETER :: PI = 4.D0*ATAN(1.D0)
      REAL*8 :: L,R1,R2,THETA

      IF (EVFEM_MAP.EQ.0) THEN
        ! Identity Map
        X=XHAT
        Y=YHAT
        Z=ZHAT
      ELSEIF (EVFEM_MAP.EQ.1) THEN
        ! Interesting Reservoir Shape - 3D
        X=XHAT+3.D0*COS(0.01D0*PI*YHAT+4.D0)*COS(0.008D0*PI*ZHAT+8.D0)
        Y=YHAT-4.D0*COS(0.3D0*PI*XHAT+2.D0)*SIN(0.01D0*PI*YHAT+5.D0)*
     &       COS(0.008D0*PI*ZHAT+8.D0)
        Z=ZHAT+5.D0*COS(0.3D0*PI*XHAT+3.D0)*COS(0.01D0*PI*YHAT+6.D0)*
     &       SIN(0.008D0*PI*ZHAT+9.D0)
      ELSEIF (EVFEM_MAP.EQ.2) THEN
        ! "S" shape curve
        R1=5.D0
        R2=10.D0
        L=PI/2*(R1+R2)
        IF (XHAT.LT.L/2) THEN
          THETA=PI/2*(1-2*XHAT/L)
          X=YHAT*COS(THETA)
          Y=YHAT*SIN(THETA)
          Z=ZHAT
        ELSE
          THETA=PI/2*(1+2*XHAT/L)
          X=(R1+R2-YHAT)*COS(THETA)+R1+R2
          Y=(R1+R2-YHAT)*SIN(THETA)
          Z=ZHAT
        ENDIF
      ELSEIF (EVFEM_MAP.EQ.3) THEN
        ! Convergence case in 2D for smooth grid distortion
        X=XHAT+0.03D0*COS(2.3D0*PI*XHAT)*COS(2.3D0*PI*YHAT)
        Y=YHAT-0.04D0*COS(2.3D0*PI*XHAT)*COS(2.3D0*PI*YHAT)
        Z=ZHAT
      ELSE
        STOP 'Unknown EVFEM_MAP'
      ENDIF

      END SUBROUTINE MAP

!----------------------------------------------------------------------
! bag8
!   Calculate area of star shaped N-gon in 3D by summing up area
!   of N triangles.
!
      SUBROUTINE AREA_NGON(N,X,Y,Z,A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(IN) :: X(N),Y(N),Z(N)
      REAL*8, INTENT(OUT) :: A
      REAL*8 :: X1(3),X2(3),X3(3)
      INTEGER :: I
      REAL*8, EXTERNAL :: TRI_AREA

      X1=[SUM(X),SUM(Y),SUM(Z)]/N   ! centroid, must be in convex hull

      A=0.D0
      DO I=1,N
        X2=[X(I),Y(I),Z(I)]
        IF (I.LT.N) THEN
          X3=[X(I+1),Y(I+1),Z(I+1)]
        ELSE
          X3=[X(1),Y(1),Z(1)]
        ENDIF
        A=A+TRI_AREA(X1,X2,X3)
      ENDDO

      END SUBROUTINE AREA_NGON

!----------------------------------------------------------------------
! bag8
!   Calculate normal to face and half mesh width on general distorted
!   hexahedra.
!
      SUBROUTINE COMP_NORMAL_AND_DX(XA,KDA,NMLA,DA)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: XA(3,8)   ! Coordinates of hexahedra
      INTEGER, INTENT(IN) :: KDA      ! Face number
      REAL*8, INTENT(OUT) :: NMLA(3)  ! Positive face unit normal
      REAL*8, INTENT(OUT) :: DA       ! Half mesh width
      REAL*8 :: EC(3),                ! Element centroid
     &          FC(3)                 ! Face centroid

      EC(1)=0.125D0*SUM(XA(1,:))
      EC(2)=0.125D0*SUM(XA(2,:))
      EC(3)=0.125D0*SUM(XA(3,:))

      IF (KDA.EQ.4) THEN       ! Low x face
        FC(1)=0.25D0*(XA(1,1)+XA(1,5)+XA(1,8)+XA(1,4))
        FC(2)=0.25D0*(XA(2,1)+XA(2,5)+XA(2,8)+XA(2,4))
        FC(3)=0.25D0*(XA(3,1)+XA(3,5)+XA(3,8)+XA(3,4))
      ELSEIF (KDA.EQ.5) THEN   ! Low y face
        FC(1)=0.25D0*(XA(1,1)+XA(1,2)+XA(1,6)+XA(1,5))
        FC(2)=0.25D0*(XA(2,1)+XA(2,2)+XA(2,6)+XA(2,5))
        FC(3)=0.25D0*(XA(3,1)+XA(3,2)+XA(3,6)+XA(3,5))
      ELSEIF (KDA.EQ.6) THEN   ! Low z face
        FC(1)=0.25D0*(XA(1,1)+XA(1,2)+XA(1,3)+XA(1,4))
        FC(2)=0.25D0*(XA(2,1)+XA(2,2)+XA(2,3)+XA(2,4))
        FC(3)=0.25D0*(XA(3,1)+XA(3,2)+XA(3,3)+XA(3,4))
      ELSEIF (KDA.EQ.1) THEN   ! High x face
        FC(1)=0.25D0*(XA(1,2)+XA(1,6)+XA(1,7)+XA(1,3))
        FC(2)=0.25D0*(XA(2,2)+XA(2,6)+XA(2,7)+XA(2,3))
        FC(3)=0.25D0*(XA(3,2)+XA(3,6)+XA(3,7)+XA(3,3))
      ELSEIF (KDA.EQ.2) THEN   ! High y face
        FC(1)=0.25D0*(XA(1,4)+XA(1,3)+XA(1,7)+XA(1,8))
        FC(2)=0.25D0*(XA(2,4)+XA(2,3)+XA(2,7)+XA(2,8))
        FC(3)=0.25D0*(XA(3,4)+XA(3,3)+XA(3,7)+XA(3,8))
      ELSEIF (KDA.EQ.3) THEN   ! High z face
        FC(1)=0.25D0*(XA(1,5)+XA(1,6)+XA(1,7)+XA(1,8))
        FC(2)=0.25D0*(XA(2,5)+XA(2,6)+XA(2,7)+XA(2,8))
        FC(3)=0.25D0*(XA(3,5)+XA(3,6)+XA(3,7)+XA(3,8))
      ENDIF

      NMLA=FC-EC
      DA=SQRT(NMLA(1)**2+NMLA(2)**2+NMLA(3)**2)
      NMLA=ABS(NMLA/DA)
      DA=2.D0*DA

      END SUBROUTINE COMP_NORMAL_AND_DX

