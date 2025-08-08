C  BLKSIN.DF - READ MULTIBLOCK DATA AND DO INITIAL PROCESSING

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE SBLKIN (NERR)
C  SUBROUTINE SETINF (IDA,JDA,KDA,IL1A,IL2A,JL1VA,JL2VA,KL1A,KL2A,
C                     KEYA,NA,NB,NIF)
C  SUBROUTINE DUALINIT (NERR)
C  SUBROUTINE IFDAT(NA,NB,IA,JA,KA,IB,JB,KB,LDIRS,A)
C  SUBROUTINE BLKFACES(NERR)
C  SUBROUTINE GETFACEINF(NBLKS,NINTFS,NERR)
C  SUBROUTINE GETFACEEL(BLK,N1,N2,X1,X2,NN1,NN2,XN1,XN2,NBLK,FACE)
C  SUBROUTINE GETINITPROJ(NBLKS,NINTFS,NUMF,FGID,FOID,NERR)
C  SUBROUTINE GETPROJMAT(NBLKS,NINTFS,NUMF,FGID,FOID,NERR)
C  SUBROUTINE GETFACEINITPRJ(BLK,N1,N2,X1,X2,NN1,NN2,XN1,XN2,NBLK
C                           ,FACE)
C  SUBROUTINE GETFACEPRJMAT(BLK,N1,N2,X1,X2,NN1,NN2,XN1,XN2,NBLK
C                          ,FACE)
C  SUBROUTINE ITRNSDPTHRCK(NERR)
C  SUBROUTINE INTFDEPTH (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
C                         KL1,KL2,KEYOUT,NBLK,DEPTH)
C  SUBROUTINE INTFROCK (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
C                       KL1,KL2,KEYOUT,NBLK,KROCK)
C  SUBROUTINE INTFTRANS (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
C                        KL1,KL2,KEYOUT,NBLK,TX,TY,TZ)
C  SUBROUTINE DUALQUIT()

C  CODE HISTORY:

C  JOHN WHEELER      2/15/99    ALPHA CODE
C  JOHN WHEELER      3/14/99    MEMORY MANAGEMENT FOR INTERFACE BUFFERS
C  SUNIL THOMAS      8/8/07     CORRECTED LATENT BUG THAT ASSIGNS NERRI
C                               TO UNINITIALIZED LOCAL NERR IN SETINF
C  SUNIL THOMAS      6/18/08    INCLUDED DATA ARRAYS FOR EASY HANDLING
C                  - 9/31/09    OF DIFFUSION-DISPERSION TERM IN REACTIVE
C                               TRANSPORT AS WELL AS THERMAL CONDUCTION
C                               AND POSSIBLY OTHER APPS. OTHER BUG FIXES
C                               NOTE: MAY NEED TO ACCOUNT FOR UNKMAP FOR
C                               TRMODEL INDICES IN MULTI-MODEL CASES
C  BEN GANIS         12/2/15    LOCAL FLUX IMPLEMENTATION FOR MULTIBLOCK
C                               WITH GENERAL HEXAHEDRA

C*********************************************************************
      SUBROUTINE SBLKIN (NERR)
C*********************************************************************

C  Reads multiblock data

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************
      USE dualmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'blkary.h'
      INCLUDE 'sblkc.h'
      INTEGER I,K,NDUM,NERR
      LOGICAL, SAVE :: ONCEONLY = .TRUE.

      NFACES=0

C INPUT AND VALIDATE INTERFACE BLOCKS

      IF (LEVELC) THEN
         WRITE (NFOUT,*)
         TITU='******'
         CALL PRTTIT(TITU)
         TITU='BLOCK INTERFACE DATA'
         CALL PRTTIT(TITU)
      ENDIF

! bag8 - added ONCEONLY for adaptivity case
      IF (ONCEONLY) THEN
      DO 1 I=1,11
      NIBLKS(1,I)=-99
    1 NIBLKS(2,I)=-99
      CALL GETVAL('FACEBLOCKS ',NIBLKS,'I4',2,11,0,0,NDUM,NERR)

!      IF (DEALII) THEN
!      NDUM=2*(NUMBLK-1)
!      DO I=1,NUMBLK-1
!      NIBLKS(1,I)=I
!      NIBLKS(2,I)=I+1
!      ENDDO
!      ENDIF

      IF (NDUM.EQ.0) THEN
         IF (LEVELC) WRITE(NFOUT,2)
    2    FORMAT(' NO FAULT BLOCK INTERFACES')
         RETURN
      ENDIF
      ENDIF

      IF (LEVELC) WRITE (NFOUT,3) NDUM/2
    3 FORMAT(/I4,' INTERFACES BETWEEN FAULT BLOCKS AS FOLLOWS:')
      DO 4 I=1,11
      IF (NIBLKS(1,I).EQ.-99.AND.NIBLKS(2,I).EQ.-99) GO TO 8
      NFACES=I
      IF (LEVELC) WRITE (NFOUT,5) I,NIBLKS(1,I),NIBLKS(2,I)
    5 FORMAT(10X,'FACE',I4,', BLOCKS',I4,' AND',I4)
      IF (NIBLKS(1,I).LT.1.OR.NIBLKS(1,I).GT.NUMBLK) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,6) NIBLKS(1,I)
    6    FORMAT (/' ERROR # 501; INVALID FAULT BLOCK NUMBER',I5,
     &      ' IN INTERFACE DATA')
      ENDIF
      IF (NIBLKS(2,I).LT.1.OR.NIBLKS(2,I).GT.NUMBLK) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,6) NIBLKS(2,I)
      ENDIF
      IF (NIBLKS(1,I).EQ.NIBLKS(2,I)) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,7) NIBLKS(1,I)
    7    FORMAT (/' ERROR # 502; FAULT BLOCK NUMBER',I5,
     &      ' REPEATED IN A PAIR')
      ENDIF
    4 CONTINUE

C  INPUT AND VALIDATE INTERFACE COMMON POINTS

    8 IF (LEVELC) WRITE (NFOUT,22)
   22 FORMAT(/' COMMON TIE POINTS FOR INTERFACES:'/
     & ' FACE      Xa          Ya          Za         ',
     & 'Xb          Yb          Zb')

      IF (ONCEONLY) THEN
      DO 9 I=1,6
      DO 9 K=1,NFACES
    9 FACECP(I,K)=-98.763D0
      CALL GETVAL('FACEXYZ ',FACECP,'R8',6,NFACES,0,0,NDUM,NERR)
!      IF (DEALII) FACECP=0.D0
      DO 10 K=1,NFACES
      DO 11 I=1,6
      IF (ABS(FACECP(I,K)+98.763D0).LT..001D0) GO TO 14
   11 CONTINUE
      GO TO 12
   14 NERR=NERR+1
      IF (LEVELC) WRITE (NFOUT,15) K
   15 FORMAT (/' ERROR # 503; INSUFICIENT DATA INPUT TO DEFINE',
     &   ' INTERFACE',I5,' COMMON POINT(S)')
      GO TO 16
   12 IF (LEVELC) THEN
         WRITE (NFOUT,17) K,(FACECP(I,K),I=1,6)
   17    FORMAT(I5,6G12.5)
      ENDIF
   16 DO 18 I=1,3
      IF (ABS(DOWN(I,NIBLKS(1,K))-DOWN(I,NIBLKS(2,K))).LT..0001)GO TO 10
   18 CONTINUE
      NERR=NERR+1
      IF (LEVELC) WRITE (NFOUT,19) NIBLKS(1,K),NIBLKS(2,K)
   19    FORMAT (/' ERROR # 504; NO PARALLEL GRID FACES BETWEEN',
     &      ' FAULT BLOCKS',I4,' AND',I4)
   10 CONTINUE
      ELSE
        DO K=1,NFACES
          IF (LEVELC) WRITE (NFOUT,17) K,(FACECP(I,K),I=1,6)
        ENDDO
      ENDIF

! bag8 - added ONCEONLY for adaptivity case
      ONCEONLY=.FALSE.

Cbag8 - read in hard coded mapping when EVFEM_HEX > 1
C    (See: subroutine MAP in infcomm.df)
        EVFEM_MAP=0
C      IF (EVFEM_HEX.GT.1) THEN
C        CALL GETVAL('EVFEM_MAP ',EVFEM_MAP,'I4',0,0,0,0,NDUM,NERR)
C        WRITE(NFOUT,*)'EVFEM_MAP=',EVFEM_MAP
C      ENDIF

! bag8 - moved to ipars.df
!      CALL DUALINIT(NERR)

      IF (NERR.NE.0) STOP 'Problem in SBLKIN'

      END
C*********************************************************************
      SUBROUTINE DUALINIT (NERR)
C*********************************************************************

C  Sets up the fault block interface system

C  NERR = Error number steped by 1 on error (input & output, INTEGER)

C*********************************************************************
      USE dualmod
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'blkary.h'
      INCLUDE 'unitsex.h'
      INCLUDE 'sblkc.h'
      INCLUDE 'mpif.h'

      INTEGER I,J,K,L,M,N,KERR,NERR
      INTEGER NFESND
      REAL*8 RNEIFS(11)
      LOGICAL :: ONCEONLY = .TRUE.     ! bag8

C  INITIALIZE POINTERS AND COUNTERS

      CALL DEALLOC_DUALMOD()
      NIFNOD=0
      NIFACE=0
      NORMFACES=0

      NKDIRS=6
      NFIEBS=1
      NFICGES=1
      NFESND=1
      DO 1 I=1,NUMBLK
      IIEBS(I)=0
    1 NIEBS(I)=0
      DO 9 I=1,11
      NEIFS(I)=0
    9 AREAT(I)=0.D0

C  DETERMINE INTERFACE ELEMENTS
C  GROUP DATA FOR EACH FAULT BLOCK TOGETHER

      NERRI=NERR
      DO 2 M=1,NUMBLK
      DO 2 N=1,NFACES
      IF (M.EQ.NIBLKS(1,N)) THEN
         CALL CALLFACE(N,NIBLKS(1,N),NIBLKS(2,N))
         IF (NERRI.GT.0) GO TO 4
      ENDIF
      IF (M.EQ.NIBLKS(2,N)) THEN
         CALL CALLFACE(N,NIBLKS(2,N),NIBLKS(1,N))
         IF (NERRI.GT.0) GO TO 4
      ENDIF
    2 CONTINUE
    4 NERR=NERRI
      IF (NERR.NE.0) RETURN

      CALL SORTIF()

C  OUTPUT SOME DATA

      CALL SUMIT(NFACES,AREAT)
      DO 15 I=1,NFACES
   15 RNEIFS(I)=NEIFS(I)+.001
      CALL SUMIT(NFACES,RNEIFS)
      DO 16 I=1,NFACES
   16 NEIFS(I)=RNEIFS(I)

      IF (LEVELC) THEN
         WRITE (NFOUT,6)
    6    FORMAT(/' FACE       AREA   INTERACTIONS')
         DO 7 I=1,NFACES
    7    WRITE(NFOUT,8) I,AREAT(I),NEIFS(I)
    8    FORMAT(I5,F12.2,I10)
         WRITE(NFOUT,*)'NIFNOD = ',NIFNOD
         WRITE(NFOUT,*)'NIFACE = ',NIFACE
      ENDIF

! bag8 debug
!      WRITE(*,*)'CIFNOD ='
!      DO I=1,NIFNOD
!        WRITE(*,'(F11.4,A,F11.4,A,F11.4)')
!     &    CIFNOD(I,1),',',CIFNOD(I,2),',',CIFNOD(I,3)
!      ENDDO
!      STOP 57
! bag8 debug

C  bag8 : allocate dual arrays for MFMFE
      IF (EVFEM_HEX.EQ.6) THEN
        CALL ALLOC_DUALMOD()
        CALL CALLWORK(FILL_DNOD2FAC,0)
      ENDIF

C  SETUP INTERFACE COMMUNICATIONS

      CALL IFCOMM (NERR)
      IF (NERR.NE.0) RETURN

C  MAP PRIMARY VARIABLES BETWEEN PHYSICAL MODELS

      DO 11 I=1,19
      DO 11 J=1,19
      DO 11 K=1,MXNUMEQS
      DO 11 L=1,MXNUMEQS
   11 UNKMAP(L,K,J,I)=0.D0

C  SAME PHYSICAL MODEL IN BOTH BLOCKS

      DO 12 I=1,19
      DO 12 K=1,MODEQS(I)
   12 UNKMAP(K,K,I,I)=1.D0

C  HYDROLOGY IMPLICIT AND SINGLE PHASE IMPLICIT

      UNKMAP(1,1,13,5)=1.D0
      UNKMAP(1,1,5,13)=1.D0

C  HYDROLOGY IMPLICIT AND BLACK 0IL IMPLICIT

      UNKMAP(1,1,2,5)=1.D0
      UNKMAP(1,1,5,2)=1.D0
      UNKMAP(2,2,2,5)=STDENO
      UNKMAP(2,2,5,2)=1.D0/STDENO

C  BLACK OIL IMPLICIT AND SINGLE PHASE IMPLICIT

      UNKMAP(1,1,13,2)=1.D0
      UNKMAP(1,1,2,13)=1.D0

C  ALLOCATE INTERFACE BUFFERS AND POINTER TO THE FIRST DIMENSION

      IF (ONCEONLY) THEN
        CALL ALCIBUF (34,NFESR,N_BUFIF,KERR)
        IF (KERR.EQ.0) CALL PNTVAR (NFESR,N_BUFDIM,KERR)
        ONCEONLY = .FALSE.
      ELSE
        CALL REALCIBUF (34,NFESR,N_BUFIF,KERR)
        IF (KERR.EQ.0) CALL REPNTVAR (NFESR,N_BUFDIM,KERR)
      ENDIF
      IF (KERR.GT.0) THEN
         NERR=NERR+1
         IF (LEVELC) WRITE (NFOUT,14) KERR
   14    FORMAT (/' ERROR #',I4,'; MEMORY MANAGEMENT ERROR ALLOCATING',
     &   ' INTERFACE BUFFERS')
         RETURN
      ENDIF

C  COMPLETE DEFINITION OF INTERFACE TRANSMISSABILITY

      CALL IFTRAN()

C  SET UP COMMUNICATION FOR MULTIGRID LINEAR SOLVER

C      CALL MDUALS(NERR)

! bag8
!      WRITE(*,*)'End of DUALINIT'
! 100  GOTO 100
!      CALL MPI_ABORT(MPI_COMM_WORLD,KERR,NERR)
!      STOP 1

      END

C*********************************************************************
      SUBROUTINE SETINF(
     &    IDA,JDA,KDA,IL1A,IL2A,JL1VA,JL2VA,KL1A,KL2A,KEYA,
     &    IDB,JDB,KDB,IL1B,IL2B,JL1VB,JL2VB,KL1B,KL2B,KEYB,
     &    NA,NB,NIF)
C*********************************************************************

C  Sets up the fault block interface system.  This is NOT a work routine.
C  Called only by callface().

C*********************************************************************
      USE dualmod
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'sblkc.h'

      INTEGER IDA,JDA,KDA,IL1A,IL2A,KL1A,KL2A
      INTEGER JL1VA(KDA),JL2VA(KDA),KEYA(IDA,JDA,KDA)
      INTEGER IDB,JDB,KDB,IL1B,IL2B,KL1B,KL2B
      INTEGER JL1VB(KDB),JL2VB(KDB),KEYB(IDB,JDB,KDB)
      INTEGER NA,NB,NIF
      INTEGER I,K,N,NP,IAA,IBB,IOFFA,JOFFA,KOFFA,NXB,NYB,NZB,MERR,
     & JL1A,JL2A,IA,JA,KA,LA,IB,JB,KB,LB,LDIRS,NERR,MPB,JMNA,JMXA,
     & INOD,IOFFB,JOFFB,KOFFB,IBL,JBL,KBL,NXA,NYA,NZA

      REAL*8 DUM1,DUM2,XA,YA,ZA,XB,YB,ZB,AX(4),AY(4),BX(4),BY(4),
     & PX(9),PY(9),XA1,YA1,ZA1,XA2,YA2,ZA2,XB1,YB1,ZB1,XB2,YB2,ZB2,TOL,
     & CT11,CT12,CT13,CT21,CT22,CT23,CT31,CT32,CT33,CO1,CO2,CO3,
     & XB1B,XB2B,YB1B,YB2B,ZB1B,ZB2B,XBMN,XBMX,YBMN,YBMX,AREAB,A
      LOGICAL XYF,XZF,YZF

! bag8
      INTEGER M,M2(4)
      REAL*8 D,X2,Y2,Z2,PXMAP(9),PYMAP(9),PZMAP(9),
     &       AREAFACE,AXMAP(9),AYMAP(9),AZMAP(9)
      REAL*8, PARAMETER :: DTOL = 1.D-6
      CHARACTER*6 STR

!      WRITE(*,*)REPEAT('-',72)
!      WRITE(*,*)'In SETINF, A block=',NA,', B block=',NB

C  SETUP FOR FAULT BLOCK REVERSAL

      IF (NA.EQ.NIBLKS(1,NIF)) THEN
         IAA=1
         IBB=4
      ELSE
         IAA=4
         IBB=1
      ENDIF

C  SET COINCIDING POINT COORDINATES AND INDEX OFFSETS

      XA1=FACECP(IAA,NIF)
      YA1=FACECP(IAA+1,NIF)
      ZA1=FACECP(IAA+2,NIF)
      XB1=FACECP(IBB,NIF)
      YB1=FACECP(IBB+1,NIF)
      ZB1=FACECP(IBB+2,NIF)

      XA2=DOWN(1,NA)
      YA2=DOWN(2,NA)
      ZA2=DOWN(3,NA)
      XB2=DOWN(1,NB)
      YB2=DOWN(2,NB)
      ZB2=DOWN(3,NB)

      CALL BLKOFF(NA,IOFFA,JOFFA,KOFFA,MERR)
      CALL BLKDIM(NA,NXA,NYA,NZA,MERR)
      CALL BLKOFF(NB,IOFFB,JOFFB,KOFFB,MERR)
      CALL BLKDIM(NB,NXB,NYB,NZB,MERR)

C  DETERMINE PARALLEL GRID PLAINS
C  INSURE ADEQUATE GRID ALLOCATION WAS SPECIFIED

      XYF=.FALSE.
      XZF=.FALSE.
      YZF=.FALSE.
      IF (ABS(XA2-XB2).LT..0001D0.AND.1.GT.0) YZF=.TRUE.
      IF (ABS(YA2-YB2).LT..0001D0.AND.1.GT.0) XZF=.TRUE.
      IF (ABS(ZA2-ZB2).LT..0001D0.AND.1.GT.0) XYF=.TRUE.

C  DETERMINE TRANSFORM FROM BLOCK B COORDINATES TO BLOCK A COORDINATES

      CT11=1.D0
      CT12=0.D0
      CT13=0.D0
      CT21=0.D0
      CT22=1.D0
      CT23=0.D0
      CT31=0.D0
      CT32=0.D0
      CT33=1.D0
      CO1=XA1-XB1
      CO2=YA1-YB1
      CO3=ZA1-ZB1

      IF (.NOT.XYF.OR..NOT.XZF) THEN

         IF (XYF) THEN
            DUM1=XB2**2+YB2**2
            DUM2=(YA2*XB2-YB2*XA2)/DUM1
            DUM1=(XA2*XB2+YA2*YB2)/DUM1
            CT11=DUM1
            CT22=DUM1
            CT12=-DUM2
            CT21=DUM2
            CO1=XA1-DUM1*XB1+DUM2*YB1
            CO2=YA1-DUM2*XB1-DUM1*YB1
         ENDIF

         IF (XZF) THEN
            DUM1=XB2**2+ZB2**2
            DUM2=(ZA2*XB2-ZB2*XA2)/DUM1
            DUM1=(XA2*XB2+ZA2*ZB2)/DUM1
            CT11=DUM1
            CT33=DUM1
            CT13=-DUM2
            CT31=DUM2
            CO1=XA1-DUM1*XB1+DUM2*ZB1
            CO3=ZA1-DUM2*XB1-DUM1*ZB1
         ENDIF

         IF (YZF) THEN
            DUM1=YB2**2+ZB2**2
            DUM2=(ZA2*YB2-ZB2*YA2)/DUM1
            DUM1=(YA2*YB2+ZA2*ZB2)/DUM1
            CT22=DUM1
            CT33=DUM1
            CT23=-DUM2
            CT32=DUM2
            CO2=YA1-DUM1*YB1+DUM2*ZB1
            CO3=ZA1-DUM2*YB1-DUM1*ZB1
         ENDIF

      ENDIF

C  LOCATE XY FACE OVERLAPS

      IF (XYF) THEN
         TOL=.0001D0*(ABS(ZREC(KL2A+1+KOFFA,NA)-ZREC(KL1A+KOFFA,NA))
     &       +ABS(ZREC(NZB+1,NB)-ZREC(1,NB)))/(KL2A-KL1A+1+NZB)

         DO 1 KA=KL1A,KL2A+1
         ZA=ZREC(KA+KOFFA,NA)
         IF (KA.GT.KL2A) THEN
            JL1A=JL1VA(KL2A)
            JL2A=JL2VA(KL2A)
         ELSE
            JL1A=JL1VA(KA)
            JL2A=JL2VA(KA)
         ENDIF

         DO 1 KB=1,NZB+1
         ZB=ZREC(KB,NB)+CO3
         IF (ABS(ZA-ZB).GT.TOL) GO TO 1

         MPB=N0MAP(NB)+MIN(KB,NZB)*NYMAP(NB)
         DO 2 JB=1,NYB
         IF (PRCMAP(MPB+JB).LT.0) GO TO 2
         YB1B=YREC(JB,NB)
         YB2B=YREC(JB+1,NB)
         DO 7 IB=1,NXB
         XB1B=XREC(IB,NB)
         XB2B=XREC(IB+1,NB)

         BX(1)=CT11*XB1B+CT12*YB1B+CO1
         BX(2)=CT11*XB1B+CT12*YB2B+CO1
         BX(3)=CT11*XB2B+CT12*YB2B+CO1
         BX(4)=CT11*XB2B+CT12*YB1B+CO1
         BY(1)=CT21*XB1B+CT22*YB1B+CO2
         BY(2)=CT21*XB1B+CT22*YB2B+CO2
         BY(3)=CT21*XB2B+CT22*YB2B+CO2
         BY(4)=CT21*XB2B+CT22*YB1B+CO2
         XBMN=BX(1)
         XBMX=BX(1)
         YBMN=BY(1)
         YBMX=BY(1)
         DO 3 I=2,4
         IF (BX(I).LT.XBMN) XBMN=BX(I)
         IF (BX(I).GT.XBMX) XBMX=BX(I)
         IF (BY(I).LT.YBMN) YBMN=BY(I)
    3    IF (BY(I).GT.YBMX) YBMX=BY(I)
         AREAB=(XB2B-XB1B)*(YB2B-YB1B)

         DO 4 JA=JL1A,JL2A
         AY(1)=YREC(JA+JOFFA,NA)
         AY(2)=YREC(JA+JOFFA+1,NA)
         IF (AY(1).GE.YBMX.OR.AY(2).LE.YBMN) GO TO 4
         AY(3)=AY(2)
         AY(4)=AY(1)
         DO 5 IA=IL1A,IL2A

! bag8
         KBL=KB-KOFFB
         JBL=JB-JOFFB
         IBL=IB-IOFFB

         IF (KEYA(IA,JA,KA).EQ.1) THEN
            IF (KEYA(IA,JA,KA-1).NE.0) GO TO 5

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.1).OR.(JBL.LT.1).OR.(KBL.LT.2)) GO TO 5
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 5
            IF (ABS(KEYB(IBL,JBL,KBL-1)).NE.1.OR.KEYB(IBL,JBL,KBL).NE.0)
     &        GO TO 5
            ENDIF

            LDIRS=6
            LA=KA
            LB=KB-1
         ELSE
            IF (KEYA(IA,JA,KA-1).NE.1) GOTO 5
            IF (KEYA(IA,JA,KA).NE.0) GO TO 5

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.1).OR.(JBL.LT.1).OR.(KBL.LT.2)) GO TO 5
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 5
            IF (KEYB(IBL,JBL,KBL-1).NE.0.OR.ABS(KEYB(IBL,JBL,KBL)).NE.1)
     &        GO TO 5
            ENDIF

            LDIRS=3
            LA=KA-1
            LB=KB
         ENDIF
         AX(1)=XREC(IA+IOFFA,NA)
         AX(3)=XREC(IA+IOFFA+1,NA)
         IF (AX(1).GE.XBMX.OR.AX(3).LE.XBMN) GO TO 5
         AX(2)=AX(1)
         AX(4)=AX(3)

         CALL POLLY(AX,AY,BX,BY,NP,PX,PY,TOL)
         IF (NP.EQ.0) GO TO 5

         IF (EVFEM_HEX.GT.1) THEN
           Z2=0.5D0*(ZA+ZB)
           CALL POLLYSORT(NP,PX,PY)     ! bag8
           CALL NEW_IFACE()
           M2=[1,4,3,2] ! Reconcile order of AX with PX
           DO N=1,NP
             CALL MAP(PX(N),PY(N),Z2,PXMAP(N),PYMAP(N),PZMAP(N))
             CALL ADD_IFNOD(PXMAP(N),PYMAP(N),PZMAP(N),INOD)
             DFAC2NOD(NIFACE,N)=INOD
! bag8 - set DFAC2IJK
             DO M=1,NP
               D=SQRT((PX(N)-AX(M2(M)))**2+(PY(N)-AY(M2(M)))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,1,M,1)=IA
                 DFAC2IJK(NIFACE,1,M,2)=JA
                 DFAC2IJK(NIFACE,1,M,3)=KA
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,1,M,1)=DFAC2IJK(NIFACE,1,M,1)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,1,M,2)=DFAC2IJK(NIFACE,1,M,2)+1
                 ENDIF
               ENDIF

               D=SQRT((PX(N)-BX(M2(M)))**2+(PY(N)-BY(M2(M)))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,2,M,1)=IB-IOFFB
                 DFAC2IJK(NIFACE,2,M,2)=JB-JOFFB
                 DFAC2IJK(NIFACE,2,M,3)=KB-KOFFB
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,2,M,1)=DFAC2IJK(NIFACE,2,M,1)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,2,M,2)=DFAC2IJK(NIFACE,2,M,2)+1
                 ENDIF
               ENDIF
             ENDDO
           ENDDO
           CALL AREA_NGON(NP,PXMAP,PYMAP,PZMAP,A)
           DO N=1,4
             CALL MAP(AX(N),AY(N),Z2,AXMAP(N),AYMAP(N),AZMAP(N))
           ENDDO
           CALL AREA_NGON(4,AXMAP,AYMAP,AZMAP,AREAFACE)
           DAREARATIO(NIFACE)=A/AREAFACE
!           WRITE(*,'(a,i3,a,f9.2)')'DAREARATIO(',NIFACE,')=',
!     &       DAREARATIO(NIFACE)
         ELSE
           A=0.D0
           DO 6 N=3,NP
    6      A=A+.5D0*((PX(N)-PX(1))*(PY(N-1)-PY(1))
     &        -(PX(N-1)-PX(1))*(PY(N)-PY(1)))
           IF (A.LT..001D0*AREAB.AND.
     &        A.LT..001D0*(AX(3)-AX(1))*(AY(3)-AY(1))) GO TO 5
!           DO N=1,NP
!             CALL ADD_IFNOD(PX(N),PY(N),Z2,INOD)
!           ENDDO
         ENDIF

         IF (IAA.EQ.1) THEN
            AREAT(NIF)=AREAT(NIF)+A
            NEIFS(NIF)=NEIFS(NIF)+1
         ENDIF

CSGT         NERRI=NERR
         CALL IFDAT(NA,NB,IA,JA,LA,IB,JB,LB,LDIRS,A)
         NERR=NERRI
         IF (NERR.GT.0) GO TO 13

    5    CONTINUE
    4    CONTINUE
    7    CONTINUE
    2    CONTINUE
    1    CONTINUE
      ENDIF

C  LOCATE XZ FACE OVERLAPS

      IF (XZF) THEN

         JMNA=JL1VA(KL1A)
         JMXA=JL2VA(KL1A)
         DO 27 K=KL1A,KL2A
         IF (JL1VA(K).LT.JMNA) JMNA=JL1VA(K)
         IF (JL2VA(K).GT.JMXA) JMXA=JL2VA(K)
   27    CONTINUE
         TOL=.0001D0*(ABS(YREC(JMXA+1+JOFFA,NA)-YREC(JMNA+JOFFA,NA))
     &       +ABS(YREC(NYB+1,NB)-YREC(1,NB)))/(JMXA-JMNA+1+NYB)

         DO 21 JA=JMNA,JMXA+1
         YA=YREC(JA+JOFFA,NA)
         DO 21 JB=1,NYB+1
         YB=YREC(JB,NB)+CO2
         IF (ABS(YA-YB).GT.TOL) GO TO 21

         DO 22 KB=1,NZB
         MPB=N0MAP(NB)+KB*NYMAP(NB)+MIN(JB,NYB)
         IF (PRCMAP(MPB).LT.0) GO TO 22
         ZB1B=ZREC(KB,NB)
         ZB2B=ZREC(KB+1,NB)
         DO 29 IB=1,NXB
         XB1B=XREC(IB,NB)
         XB2B=XREC(IB+1,NB)

         BX(1)=CT11*XB1B+CT13*ZB1B+CO1
         BX(2)=CT11*XB1B+CT13*ZB2B+CO1
         BX(3)=CT11*XB2B+CT13*ZB2B+CO1
         BX(4)=CT11*XB2B+CT13*ZB1B+CO1
         BY(1)=CT31*XB1B+CT33*ZB1B+CO3
         BY(2)=CT31*XB1B+CT33*ZB2B+CO3
         BY(3)=CT31*XB2B+CT33*ZB2B+CO3
         BY(4)=CT31*XB2B+CT33*ZB1B+CO3
         XBMN=BX(1)
         XBMX=BX(1)
         YBMN=BY(1)
         YBMX=BY(1)
         DO 23 I=2,4
         IF (BX(I).LT.XBMN) XBMN=BX(I)
         IF (BX(I).GT.XBMX) XBMX=BX(I)
         IF (BY(I).LT.YBMN) YBMN=BY(I)
   23    IF (BY(I).GT.YBMX) YBMX=BY(I)
         AREAB=(XB2B-XB1B)*(YB2B-YB1B)

         DO 24 KA=KL1A,KL2A
         IF (JA.LT.JL1VA(KA).OR.JA.GT.JL2VA(KA)+1) GO TO 24
         AY(1)=ZREC(KA+KOFFA,NA)
         AY(2)=ZREC(KA+KOFFA+1,NA)
         IF (AY(1).GE.YBMX.OR.AY(2).LE.YBMN) GO TO 24
         AY(3)=AY(2)
         AY(4)=AY(1)
         DO 25 IA=IL1A,IL2A

! bag8
         KBL=KB-KOFFB
         JBL=JB-JOFFB
         IBL=IB-IOFFB

         IF (KEYA(IA,JA,KA).EQ.1) THEN
            IF (KEYA(IA,JA-1,KA).NE.0) GO TO 25

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.1).OR.(JBL.LT.2).OR.(KBL.LT.1)) GO TO 25
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 25
            IF (ABS(KEYB(IBL,JBL-1,KBL)).NE.1.OR.KEYB(IBL,JBL,KBL).NE.0)
     &        GO TO 25
            ENDIF

            LDIRS=5
            LA=JA
            LB=JB-1
         ELSE
            IF (KEYA(IA,JA-1,KA).NE.1) GOTO 25
            IF (KEYA(IA,JA,KA).NE.0) GO TO 25

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.1).OR.(JBL.LT.2).OR.(KBL.LT.1)) GO TO 25
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 25
            IF (KEYB(IBL,JBL-1,KBL).NE.0.OR.ABS(KEYB(IBL,JBL,KBL)).NE.1)
     &        GO TO 25
            ENDIF

            LDIRS=2
            LA=JA-1
            LB=JB
         ENDIF
         AX(1)=XREC(IA+IOFFA,NA)
         AX(3)=XREC(IA+IOFFA+1,NA)
         IF (AX(1).GE.XBMX.OR.AX(3).LE.XBMN) GO TO 25
         AX(2)=AX(1)
         AX(4)=AX(3)

         CALL POLLY(AX,AY,BX,BY,NP,PX,PY,TOL)
         IF (NP.EQ.0) GO TO 25

         IF (EVFEM_HEX.GT.1) THEN
           Y2=0.5D0*(YA+YB)
           CALL POLLYSORT(NP,PX,PY)    ! bag8
           CALL NEW_IFACE()
           M2=[1,4,3,2] ! Reconcile order of AX with PX
           DO N=1,NP
             CALL MAP(PX(N),Y2,PY(N),PXMAP(N),PYMAP(N),PZMAP(N))
             CALL ADD_IFNOD(PXMAP(N),PYMAP(N),PZMAP(N),INOD)
             DFAC2NOD(NIFACE,N)=INOD
! bag8 - set DFAC2IJK
             DO M=1,NP
               D=SQRT((PX(N)-AX(M2(M)))**2+(PY(N)-AY(M2(M)))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,1,M,1)=IA
                 DFAC2IJK(NIFACE,1,M,2)=JA
                 DFAC2IJK(NIFACE,1,M,3)=KA
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,1,M,1)=DFAC2IJK(NIFACE,1,M,1)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,1,M,3)=DFAC2IJK(NIFACE,1,M,3)+1
                 ENDIF
               ENDIF
               D=SQRT((PX(N)-BX(M2(M)))**2+(PY(N)-BY(M2(M)))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,2,M,1)=IB-IOFFB
                 DFAC2IJK(NIFACE,2,M,2)=JB-JOFFB
                 DFAC2IJK(NIFACE,2,M,3)=KB-KOFFB
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,2,M,1)=DFAC2IJK(NIFACE,2,M,1)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,2,M,3)=DFAC2IJK(NIFACE,2,M,3)+1
                 ENDIF
               ENDIF
             ENDDO
           ENDDO
           CALL AREA_NGON(NP,PXMAP,PYMAP,PZMAP,A)
           DO N=1,4
             CALL MAP(AX(N),Y2,AY(N),AXMAP(N),AYMAP(N),AZMAP(N))
           ENDDO
           CALL AREA_NGON(4,AXMAP,AYMAP,AZMAP,AREAFACE)
           DAREARATIO(NIFACE)=A/AREAFACE
!           WRITE(*,'(a,i3,a,f9.2)')'DAREARATIO(',NIFACE,')=',
!     &       DAREARATIO(NIFACE)
         ELSE
           A=0.D0
           DO 26 N=3,NP
   26      A=A+.5D0*((PX(N)-PX(1))*(PY(N-1)-PY(1))
     &        -(PX(N-1)-PX(1))*(PY(N)-PY(1)))
           IF (A.LT..001D0*AREAB.AND.
     &        A.LT..001D0*(AX(3)-AX(1))*(AY(3)-AY(1))) GO TO 25
!           DO N=1,NP
!             CALL ADD_IFNOD(PX(N),Y2,PY(N),INOD)
!           ENDDO
         ENDIF

         IF (IAA.EQ.1) THEN
            AREAT(NIF)=AREAT(NIF)+A
            NEIFS(NIF)=NEIFS(NIF)+1
         ENDIF

CSGT         NERRI=NERR
         CALL IFDAT(NA,NB,IA,LA,KA,IB,LB,KB,LDIRS,A)
         NERR=NERRI
         IF (NERR.GT.0) GO TO 13

   25    CONTINUE
   24    CONTINUE
   29    CONTINUE
   22    CONTINUE
   21    CONTINUE
      ENDIF

C  LOCATE YZ FACE OVERLAPS

      IF (YZF) THEN
         TOL=.0001D0*(ABS(XREC(IL2A+1+IOFFA,NA)-XREC(IL1A+IOFFA,NA))
     &       +ABS(XREC(NXB+1,NB)-XREC(1,NB)))/(IL2A-IL1A+1+NXB)

         DO 41 IA=IL1A,IL2A+1
         XA=XREC(IA+IOFFA,NA)
         DO 41 IB=1,NXB+1
         XB=XREC(IB,NB)+CO1
         IF (ABS(XA-XB).GT.TOL) GO TO 41

         DO 42 KB=1,NZB
         ZB1B=ZREC(KB,NB)
         ZB2B=ZREC(KB+1,NB)
         MPB=N0MAP(NB)+MIN(KB,NZB)*NYMAP(NB)
         DO 42 JB=1,NYB
         IF (PRCMAP(MPB+JB).LT.0) GO TO 42

         YB1B=YREC(JB,NB)
         YB2B=YREC(JB+1,NB)

         BX(1)=CT33*ZB1B+CT32*YB1B+CO3
         BX(2)=CT33*ZB1B+CT32*YB2B+CO3
         BX(3)=CT33*ZB2B+CT32*YB2B+CO3
         BX(4)=CT33*ZB2B+CT32*YB1B+CO3
         BY(1)=CT23*ZB1B+CT22*YB1B+CO2
         BY(2)=CT23*ZB1B+CT22*YB2B+CO2
         BY(3)=CT23*ZB2B+CT22*YB2B+CO2
         BY(4)=CT23*ZB2B+CT22*YB1B+CO2
         XBMN=BX(1)
         XBMX=BX(1)
         YBMN=BY(1)
         YBMX=BY(1)
         DO 43 I=2,4
         IF (BX(I).LT.XBMN) XBMN=BX(I)
         IF (BX(I).GT.XBMX) XBMX=BX(I)
         IF (BY(I).LT.YBMN) YBMN=BY(I)
   43    IF (BY(I).GT.YBMX) YBMX=BY(I)
         AREAB=(XB2B-XB1B)*(YB2B-YB1B)

         DO 44 KA=KL1A,KL2A
         AX(1)=ZREC(KA+KOFFA,NA)
         AX(3)=ZREC(KA+KOFFA+1,NA)
         IF (AX(1).GE.XBMX.OR.AX(3).LE.XBMN) GO TO 44
         AX(2)=AX(1)
         AX(4)=AX(3)
         IF (KA.GT.KL2A) THEN
            JL1A=JL1VA(KL2A)
            JL2A=JL2VA(KL2A)
         ELSE
            JL1A=JL1VA(KA)
            JL2A=JL2VA(KA)
         ENDIF
         DO 45 JA=JL1A,JL2A
         AY(1)=YREC(JA+JOFFA,NA)
         AY(2)=YREC(JA+JOFFA+1,NA)
         IF (AY(1).GE.YBMX.OR.AY(2).LE.YBMN) GO TO 45
         AY(3)=AY(2)
         AY(4)=AY(1)

! bag8
         KBL=KB-KOFFB
         JBL=JB-JOFFB
         IBL=IB-IOFFB

         IF (KEYA(IA,JA,KA).EQ.1) THEN
            IF (KEYA(IA-1,JA,KA).NE.0) GO TO 45

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.2).OR.(JBL.LT.1).OR.(KBL.LT.1)) GO TO 45
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 45
            IF (ABS(KEYB(IBL-1,JBL,KBL)).NE.1.OR.KEYB(IBL,JBL,KBL).NE.0)
     &        GO TO 45
            ENDIF

            LDIRS=4
            LA=IA
            LB=IB-1
         ELSE
            IF (KEYA(IA-1,JA,KA).NE.1) GOTO 45
            IF (KEYA(IA,JA,KA).NE.0) GO TO 45

! bag8
            IF (EV_PRCBLK.GT.0) THEN
            IF ((IBL.LT.2).OR.(JBL.LT.1).OR.(KBL.LT.1)) GO TO 45
            IF ((IBL.GT.IDB).OR.(JBL.GT.JDB).OR.(KBL.GT.KDB)) GO TO 45
            IF (KEYB(IBL-1,JBL,KBL).NE.0.OR.ABS(KEYB(IBL,JBL,KBL)).NE.1)
     &        GO TO 45
            ENDIF

            LDIRS=1
            LA=IA-1
            LB=IB
         ENDIF

         CALL POLLY(AX,AY,BX,BY,NP,PX,PY,TOL)
         IF (NP.EQ.0) GO TO 45

         IF (EVFEM_HEX.GT.1) THEN
           X2=0.5D0*(XA+XB)
           CALL POLLYSORT(NP,PY,PX)  ! bag8 - note the order is PY,PX!
           CALL NEW_IFACE()
           DO N=1,NP
             CALL MAP(X2,PY(N),PX(N),PXMAP(N),PYMAP(N),PZMAP(N))
             CALL ADD_IFNOD(PXMAP(N),PYMAP(N),PZMAP(N),INOD)
             DFAC2NOD(NIFACE,N)=INOD
! bag8 - set DFAC2IJK
             DO M=1,NP
               D=SQRT((PX(N)-AX(M))**2+(PY(N)-AY(M))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,1,M,1)=IA
                 DFAC2IJK(NIFACE,1,M,2)=JA
                 DFAC2IJK(NIFACE,1,M,3)=KA
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,1,M,2)=DFAC2IJK(NIFACE,1,M,2)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,1,M,3)=DFAC2IJK(NIFACE,1,M,3)+1
                 ENDIF
               ENDIF
               D=SQRT((PX(N)-BX(M))**2+(PY(N)-BY(M))**2)
               IF (D.LT.DTOL) THEN
                 DFAC2IJK(NIFACE,2,M,1)=IB-IOFFB
                 DFAC2IJK(NIFACE,2,M,2)=JB-JOFFB
                 DFAC2IJK(NIFACE,2,M,3)=KB-KOFFB
                 IF ((M.EQ.2).OR.(M.EQ.3)) THEN
                   DFAC2IJK(NIFACE,2,M,2)=DFAC2IJK(NIFACE,2,M,2)+1
                 ENDIF
                 IF ((M.EQ.3).OR.(M.EQ.4)) THEN
                   DFAC2IJK(NIFACE,2,M,3)=DFAC2IJK(NIFACE,2,M,3)+1
                 ENDIF
               ENDIF
             ENDDO
           ENDDO
           CALL AREA_NGON(NP,PXMAP,PYMAP,PZMAP,A)
           DO N=1,4
             CALL MAP(X2,AY(N),AX(N),AXMAP(N),AYMAP(N),AZMAP(N))
           ENDDO
           CALL AREA_NGON(4,AXMAP,AYMAP,AZMAP,AREAFACE)
           DAREARATIO(NIFACE)=A/AREAFACE
!           WRITE(*,'(a,i3,a,f9.2)')'DAREARATIO(',NIFACE,')=',
!     &       DAREARATIO(NIFACE)
         ELSE
           A=0.D0
           DO 46 N=3,NP
   46      A=A+.5D0*((PX(N)-PX(1))*(PY(N-1)-PY(1))
     &        -(PX(N-1)-PX(1))*(PY(N)-PY(1)))
           IF (A.LT..001D0*AREAB.AND.
     &        A.LT..001D0*(AX(3)-AX(1))*(AY(3)-AY(1))) GO TO 45
!           DO N=1,NP
!             CALL ADD_IFNOD(X2,PY(N),PX(N),INOD)
!           ENDDO
         ENDIF

         IF (IAA.EQ.1) THEN
            AREAT(NIF)=AREAT(NIF)+A
            NEIFS(NIF)=NEIFS(NIF)+1
         ENDIF

CSGT         NERRI=NERR
         CALL IFDAT(NA,NB,LA,JA,KA,LB,JB,KB,LDIRS,A)
         NERR=NERRI
         IF (NERR.GT.0) GO TO 13

   45    CONTINUE
   44    CONTINUE
   42    CONTINUE
   41    CONTINUE
      ENDIF

   13 CONTINUE
      END

!*********************************************************************
! bag8 - Sort the result from POLLY
!*********************************************************************
      SUBROUTINE POLLYSORT(NP,PX,PY)
!*********************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NP
      REAL*8, INTENT(INOUT) :: PX(NP), PY(NP)
      REAL*8 :: PXIN(4), PYIN(4)
      REAL*8 :: XLO,XHI,YLO,YHI,D
      REAL*8, PARAMETER :: TOL = 1.D-6
      INTEGER N

      IF (NP.NE.4) THEN
        WRITE(*,*)'NP=',NP
        STOP 'In POLLYSORT, NP not equal 4'
      ENDIF

      PXIN=PX
      PYIN=PY
      XLO=MINVAL(PX)
      XHI=MAXVAL(PX)
      YLO=MINVAL(PY)
      YHI=MAXVAL(PY)

      DO N=1,NP
        D=SQRT((PXIN(N)-XLO)**2+(PYIN(N)-YLO)**2)
        IF (D.LT.TOL) THEN
          PX(1)=PXIN(N)
          PY(1)=PYIN(N)
          GOTO 1
        ENDIF
      ENDDO
      STOP 'In CALC_NSORT, could not find node 1'
 1    CONTINUE

      DO N=1,NP
        D=SQRT((PXIN(N)-XHI)**2+(PYIN(N)-YLO)**2)
        IF (D.LT.TOL) THEN
          PX(2)=PXIN(N)
          PY(2)=PYIN(N)
          GOTO 2
        ENDIF
      ENDDO
      STOP 'In CALC_NSORT, could not find node 2'

 2    CONTINUE
      DO N=1,NP
        D=SQRT((PXIN(N)-XHI)**2+(PYIN(N)-YHI)**2)
        IF (D.LT.TOL) THEN
          PX(3)=PXIN(N)
          PY(3)=PYIN(N)
          GOTO 3
        ENDIF
      ENDDO
      STOP 'In CALC_NSORT, could not find node 3'

 3    CONTINUE
      DO N=1,NP
        D=SQRT((PXIN(N)-XLO)**2+(PYIN(N)-YHI)**2)
        IF (D.LT.TOL) THEN
          PX(4)=PXIN(N)
          PY(4)=PYIN(N)
          GOTO 4
        ENDIF
      ENDDO
      STOP 'In CALC_NSORT, could not find node 4'

 4    CONTINUE

      END SUBROUTINE POLLYSORT

C*********************************************************************
      SUBROUTINE IFDAT(NA,NB,IA,JA,KA,IB,JB,KB,LDIRS,A)
C*********************************************************************

C  Records block interface data in COMMON /SBLKC/ for one element of the
C  interface

C  NA = Fault block number, block A (input, INTEGER)

C  NB = Fault block number, block B (input, INTEGER)

C  IA,JA,KA = Local coordinates of the block A element connected to the
C             interface (input, INTEGER)

C  IB,JB,KB = Global coordinates of the block B element connected to the
C             interface (input, INTEGER)

C  LDIRS    = Direction key for element B relative to element A
C             (input, INTEGER)
C           =  1 ==> X + 1
C           =  2 ==> Y + 1
C           =  3 ==> Z + 1
C           =  4 ==> X - 1
C           =  5 ==> Y - 1
C           =  6 ==> Z - 1
C           =  Additional offsets for higher order approximations.

C  A = Interface sub area, sq-ft (input, REAL*8)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'

      INCLUDE 'sblkc.h'

      INTEGER IA,JA,KA,NA,IB,JB,KB,NB,K,N,N1,N2,LDIRS
      REAL*8 A

C  REGISTER AN A ELEMENT IF NOT PREVIOUSLY REGISTERED

      IF (IIEBS(NA).EQ.0) IIEBS(NA)=NFIEBS
      N1=IIEBS(NA)
      N2=N1+NIEBS(NA)-1
      DO 1 N=N1,N2
      IF (IA.EQ.IJKS(1,N).AND.JA.EQ.IJKS(2,N).AND.KA.EQ.IJKS(3,N)) THEN
         K=N
         GO TO 2
      ENDIF
    1 CONTINUE
      K=NFIEBS
      IF (NFIEBS.GE.40000) THEN
         NERRI=NERRI+1
         IF (LEVELC) WRITE (NFOUT,6) NFIEBS-1
    6    FORMAT (/' ERROR # 505; MAX NUMBER OF INTERFACE ELEMENTS (',
     &      I7,') EXCEEDED')
         RETURN
      ENDIF
      NFIEBS=K+1
      NIEBS(NA)=NIEBS(NA)+1
      IJKS(1,K)=IA
      IJKS(2,K)=JA
      IJKS(3,K)=KA
      NCGES(K)=0
      ICGES(K)=NFICGES

C  REGISTER AN B ELEMENT
C  INITIALIZE IN ANY ORDER AND LATER SORT BY NAEB
C  NOTE THAT ICGES(K) MAY CHANGE ON SORTING

    2 NCGES(K)=NCGES(K)+1
      IF (NFICGES.GT.90000) THEN
         NERRI=NERRI+1
         IF (LEVELC) WRITE (NFOUT,7) NFICGES-1
    7    FORMAT (/' ERROR # 506; MAX NUMBER OF INTERFACE COUPLINGS (',
     &      I7,') EXCEEDED')
         RETURN
      ENDIF
      KDIRS(NFICGES)=LDIRS
      NAEB(NFICGES)=K
      TFINS(NFICGES)=A
      AREAI(NFICGES)=A
C COPY FOR ISOLATED MATRIX SOLVES (NON-MODEL DEPENDENT)
      TDDFINS0(NFICGES)=A
C
      JBLOCK(NFICGES)=NB
      IJKT(1,NFICGES)=IB
      IJKT(2,NFICGES)=JB
      IJKT(3,NFICGES)=KB
C COPY FOR ISOLATED MATRIX SOLVES (NON-MODEL DEPENDENT)
      TRDDIJKT(1,NFICGES)=IB
      TRDDIJKT(2,NFICGES)=JB
      TRDDIJKT(3,NFICGES)=KB
C
      NFICGES=NFICGES+1

      END

C*********************************************************************

      SUBROUTINE DUALQUIT()
      USE dualmod
      IMPLICIT NONE
      CALL DEALLOC_DUALMOD()
      END
