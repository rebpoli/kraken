C  INITIAL.F - INITIALIZE DATA IN GRID-ELEMENT ARRAYS - ALL MODELS

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE DEPTH1   (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,DEPTH)
C  SUBROUTINE DEPTH2   (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,DEPTH,XC,YC,ZC)
C  SUBROUTINE DEPTH1_TOP   (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,DEPTH,XC,YC,ZC)
C  SUBROUTINE TRANC1   (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,TX,TY,TZ,PX,PY,PZ)
C  SUBROUTINE SETARYI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,IARY,I4)
C  SUBROUTINE SETARYI4P (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                       KL2,KEYOUT,NBLK,IARY,I4)
C  SUBROUTINE SETARYR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,RARY,R4)
C  SUBROUTINE SETARYR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,RARY,R8)
C  SUBROUTINE SETARYR8N (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,RARY,R8,N)
C  SUBROUTINE CPYARYI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,IARY1,IARY2)
C  SUBROUTINE CPYARYR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,RARY1,RARY2)
C  SUBROUTINE CPYARYR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                      KL2,KEYOUT,NBLK,RARY1,RARY2)
C  SUBROUTINE CPYARYR8N(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
C                       KL2,KEYOUT,NBLK,RARY1,RARY2,N)
C  SUBROUTINE COUNTE   (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                      KEYOUT,NBLK)

C  CODE HISTORY:

C  JOHN WHEELER     1/17/95     ALPHA CODE
C  JOHN WHEELER     3/28/00     COUNT GRID ELEMENTS
C  RICK DEAN        1/18/02     ADDED SETARYR8N AND CPYARYR8N
C  RICK DEAN       12/19/02     ADDED SETARYI4P
C  XIANHUI KONG    06/11/14     ADDED DEPTH1_TOP READIN TOPSURFACE

C*********************************************************************
      SUBROUTINE DEPTH1 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,DEPTH)
C*********************************************************************

C  Calculate block center depth array for orthogonal grid option.
C  This is a work routine.
C  DEPTH = 0. at x=0., y=0., z=0.

C  DEPTH(I,J,K) = Block center depth in feet for fault block NBLK
C                 (output, REAL*8)

C*********************************************************************
C      INCLUDE 'msjunk.h'

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 DEPTH(IDIM,JDIM,KDIM)

      CALL BLKDIM(NBLK,ID,JD,KD,MERR)
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

      DO 4 K=1,KDIM
      KG=K+KOFF
      IF (KG.GT.0.AND.KG.LE.KD) THEN
         DO 5 J=1,JDIM
         JG=J+JOFF
         IF (JG.GT.0.AND.JG.LE.JD) THEN
            DO 6 I=1,IDIM
            IG=I+IOFF
            IF (IG.GT.0.AND.IG.LE.ID) DEPTH(I,J,K)=.5D0*(
     &         DOWN(1,NBLK)*(XREC(IG,NBLK)+XREC(IG+1,NBLK))+
     &         DOWN(2,NBLK)*(YREC(JG,NBLK)+YREC(JG+1,NBLK))+
     &         DOWN(3,NBLK)*(ZREC(KG,NBLK)+ZREC(KG+1,NBLK)))
    6       CONTINUE
         ENDIF
    5    CONTINUE
      ENDIF
    4 CONTINUE

      END
C*********************************************************************
      SUBROUTINE DEPTH2 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK,DEPTH,XC,YC,ZC)
C*********************************************************************

C  Calculate block center depth array for corner point grid option.
C  This is a work routine.
C  DEPTH = 0. at x=0., y=0., z=0.

C  DEPTH(I,J,K) = Block center depth in feet for fault block NBLK
C                 (output, REAL*8)

C  XC(I,J,K),YC(I,J,K),ZC(I,J,K) = Corner point locations for fault
C                 block NBLK (input, REAL*4)

C*********************************************************************

      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      REAL*8 X,Y,Z,DEPTH(IDIM,JDIM,KDIM)
      REAL*4 XC(IDIM,JDIM,KDIM),YC(IDIM,JDIM,KDIM),ZC(IDIM,JDIM,KDIM)

      DO 1 K=2,KDIM
      DO 1 J=2,JDIM
      DO 1 I=2,IDIM
      X=XC(I,J,K)
      X=0.125D0*(X+XC(I,J,K-1)+XC(I,J-1,K)+XC(I,J-1,K-1)+
     &   XC(I-1,J,K)+XC(I-1,J,K-1)+XC(I-1,J-1,K)+XC(I-1,J-1,K-1))
      Y=YC(I,J,K)
      Y=0.125D0*(Y+YC(I,J,K-1)+YC(I,J-1,K)+YC(I,J-1,K-1)+
     &   YC(I-1,J,K)+YC(I-1,J,K-1)+YC(I-1,J-1,K)+YC(I-1,J-1,K-1))
      Z=ZC(I,J,K)
      Z=0.125D0*(Z+ZC(I,J,K-1)+ZC(I,J-1,K)+ZC(I,J-1,K-1)+
     &   ZC(I-1,J,K)+ZC(I-1,J,K-1)+ZC(I-1,J-1,K)+ZC(I-1,J-1,K-1))
    1 DEPTH(I-1,J-1,K-1)=DOWN(1,NBLK)*X+DOWN(2,NBLK)*Y+DOWN(3,NBLK)*Z

      END
C*********************************************************************
      SUBROUTINE TRANC1(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &   KL1,KL2,KEYOUT,NBLK,TX,TY,TZ,PX,PY,PZ)
C*********************************************************************

C  Calculate transmissability constant array for orthogonal grid option.
C  This is a work routine.

C  TX(I,J,K) = Transmissability (OUTPUT, REAL*8)
C  TY(I,J,K)
C  TZ(I,J,K)

C  PX(I,J,K) = Permeability (INPUT AND OUTPUT, REAL*4)
C  PY(I,J,K)   (0 is put in keyed out elements)
C  PZ(I,J,K)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'layout.h'
      INCLUDE 'blkary.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1,JL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 TX(IDIM,JDIM,KDIM),TY(IDIM,JDIM,KDIM),TZ(IDIM,JDIM,KDIM)
      REAL*4 PX(IDIM,JDIM,KDIM),PY(IDIM,JDIM,KDIM),PZ(IDIM,JDIM,KDIM)
      INTEGER I,J,K,IOFF,JOFF,KOFF,MERR,ID,JD,KD
      REAL*8 DX,DY,DZ,DUM1,DUM2,CVC

C bag8 - CVF = 2 sq-ft cp / md psi day
      CVC = 2 * CONV_FACTOR

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      CALL BLKDIM(NBLK,ID,JD,KD,MERR)

C  PUT 0. PERMEABILITY IN UNUSED LOCATIONS

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
      TX(I,J,K)=0.D0
      TY(I,J,K)=0.D0
      TZ(I,J,K)=0.D0
! bag8 - do not zero out perm for adaptivity to work
!      IF (KEYOUT(I,J,K).EQ.0) THEN
!         PX(I,J,K)=0.
!         PY(I,J,K)=0.
!         PZ(I,J,K)=0.
!      ENDIF
    1 CONTINUE

C  CONVERT FROM PERMEABILITIES TO TRANSMISSABILITY COEFFICIENT CONSTANT

      DO 2 K=KL1,KL2
      DZ=DZREC(K+KOFF,NBLK)
      JL1=JL1V(K)
      JL2=JL2V(K)
      DO 2 J=JL1,JL2
      DY=DYREC(J+JOFF,NBLK)
      IF (IL1.GT.1) THEN
         DUM1=PX(IL1-1,J,K)
      ELSE
         DUM1=0.D0
      ENDIF
      DO 3 I=IL1,IL2
      DX=DXREC(I+IOFF,NBLK)
      DUM2=PX(I,J,K)
      TX(I,J,K)=0.D0
      IF (I+IOFF.GT.1) THEN
        IF (DUM1.GT.0.D0.AND.DUM2.GT.0.D0.AND.
     &   (KEYOUT(I,J,K).NE.0.AND.KEYOUT(I-1,J,K).NE.0).AND.
     &   (I+IOFF-1.GT.0)) THEN
         TX(I,J,K)=CVC*DY*DZ/(DXREC(I+IOFF-1,NBLK)/DUM1+DX/DUM2)
        ENDIF
      ENDIF
    3 DUM1=DUM2
      IF (DUM2.GT.0.D0.AND.IL2.LT.IDIM.AND.PX(IL2+1,J,K).GT.0..AND.
     & (KEYOUT(IL2+1,J,K).NE.0.AND.KEYOUT(IL2,J,K).NE.0).AND.
     & (IL2+IOFF+1.LE.ID+1)) THEN
         TX(IL2+1,J,K)=CVC*DY*DZ/(DXREC(IL2+IOFF+1,NBLK)/PX(IL2+1,J,K)
     &    +DX/DUM2)
      ELSE
         IF (IL2.LT.IDIM) TX(IL2+1,J,K)=0.D0
      ENDIF
    2 CONTINUE

      DO 4 K=KL1,KL2
      DZ=DZREC(K+KOFF,NBLK)
      JL1=JL1V(K)
      JL2=JL2V(K)
      DO 4 I=IL1,IL2
      DX=DXREC(I+IOFF,NBLK)
      IF (JL1.GT.1) THEN
         DUM1=PY(I,JL1-1,K)
      ELSE
         DUM1=0.D0
      ENDIF
      DO 5 J=JL1,JL2
      IF ((J+JOFF.LT.1).OR.(J+JOFF.GT.JD)) CYCLE
      DY=DYREC(J+JOFF,NBLK)
      DUM2=PY(I,J,K)
      IF (DUM1.GT.0.D0.AND.DUM2.GT.0.D0.AND.J.GT.1.AND.
     & (KEYOUT(I,J,K).NE.0.AND.KEYOUT(I,J-1,K).NE.0).AND.
     & (J+JOFF-1.GT.0)) THEN
         TY(I,J,K)=CVC*DX*DZ/(DYREC(J+JOFF-1,NBLK)/DUM1+DY/DUM2)
      ELSE
         TY(I,J,K)=0.D0
      ENDIF
    5 DUM1=DUM2
      IF (DUM2.GT.0.D0.AND.JL2.LT.JDIM.AND.PY(I,JL2+1,K).GT.0..AND.
     & (KEYOUT(I,JL2+1,K).NE.0.AND.KEYOUT(I,JL2,K).NE.0).AND.
     & (JL2+JOFF+1.LE.JD+1)) THEN
         TY(I,JL2+1,K)=CVC*DX*DZ/(DYREC(JL2+JOFF+1,NBLK)/PY(I,JL2+1,K)
     &    +DY/DUM2)
      ELSE
         TY(I,JL2+1,K)=0.D0
      ENDIF
    4 CONTINUE

      JL1=JDIM
      JL2=1
      DO 6 K=KL1,KL2
      IF (JL1V(K).LT.JL1) JL1=JL1V(K)
      IF (JL2V(K).GT.JL2) JL2=JL2V(K)
    6 CONTINUE

      DO 7 J=JL1,JL2
      IF ((J+JOFF.LT.1).OR.(J+JOFF.GT.JD)) CYCLE
      DY=DYREC(J+JOFF,NBLK)
      DO 7 I=IL1,IL2
      IF ((I+IOFF.LT.1).OR.(I+IOFF.GT.ID)) CYCLE
      DX=DXREC(I+IOFF,NBLK)
      IF (KL1.GT.1) THEN
         DUM1=PZ(I,J,KL1-1)
      ELSE
         DUM1=0.D0
      ENDIF
      DO 8 K=KL1,KL2
      DZ=DZREC(K+KOFF,NBLK)
      DUM2=PZ(I,J,K)
      IF (DUM1.GT.0.D0.AND.DUM2.GT.0.D0.AND.K.GT.1.AND.
     & (KEYOUT(I,J,K).NE.0.AND.KEYOUT(I,J,K-1).NE.0).AND.
     & (K+KOFF-1.GT.0)) THEN
         TZ(I,J,K)=CVC*DX*DY/(DZREC(K+KOFF-1,NBLK)/DUM1+DZ/DUM2)
      ELSE
         TZ(I,J,K)=0.D0
      ENDIF
    8 DUM1=DUM2
      IF (DUM2.GT.0.D0.AND.KL2.LT.KDIM.AND.PZ(I,J,KL2+1).GT.0..AND.
     & (KEYOUT(I,J,KL2+1).NE.0.AND.KEYOUT(I,J,KL2).NE.0).AND.
     & (KL2+KOFF+1.LE.KD+1)) THEN
         TZ(I,J,KL2+1)=CVC*DX*DY/(DZREC(KL2+KOFF+1,NBLK)/PZ(I,J,KL2+1)
     &    +DZ/DUM2)
      ELSE
         IF (KL2.LT.KDIM) TZ(I,J,KL2+1)=0.D0
      ENDIF
    7 CONTINUE

      END
C*********************************************************************
      SUBROUTINE SETARYI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,IARY,I4)
C*********************************************************************

C  Set ALL elements of an INTEGER grid-element array to a constant.
C  This is a work routine.

C  IARY(I,J,K) = Array to be initialized (output, INTEGER)

C  I4 = Value to which the array elements are to be set
C       (input, INTEGER)

C*********************************************************************

      INTEGER I4,IARY(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 IARY(I,J,K)=I4

      END
C*********************************************************************
      SUBROUTINE SETARYI4P (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                      KL2,KEYOUT,NBLK,IARY,I4)
C*********************************************************************

C  Set active elements of an INTEGER grid-element array to a constant.
C  Inactive elements set to zero.
C  This is a work routine.

C  IARY(I,J,K) = Array to be initialized (output, INTEGER)

C  I4 = Value to which the array elements are to be set
C       (input, INTEGER)

C*********************************************************************

      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM*JDIM*KDIM)
      INTEGER I4,IARY(IDIM*JDIM*KDIM)

      DO 1 K=1,IDIM*JDIM*KDIM
         IF(KEYOUT(K).EQ.1) THEN
            IARY(K)=I4
         ELSE
            IARY(K)=0
         ENDIF
    1 CONTINUE

      END
C*********************************************************************
      SUBROUTINE SETARYR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,RARY,R4)
C*********************************************************************

C  Set ALL elements of an REAL*4 grid-element array to a constant.
C  This is a work routine.

C  RARY(I,J,K) = Array to be initialized (output, REAL*4)

C  R4 = Value to which the array elements are to be set
C       (input, REAL*4)

C*********************************************************************

      REAL*4 R4,RARY(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 RARY(I,J,K)=R4

      END

C*********************************************************************
      SUBROUTINE SETARYR4N (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                      KL2,KEYOUT,NBLK,RARY,R4,N)
C*********************************************************************

C  Set ALL elements of a multidimensional REAL*4 grid-element array
C  to a constant.
C  This is a work routine.

C  RARY(I,J,K,L) = Array to be initialized (output, REAL*8)

C  R4 = Value to which the array elements are to be set
C       (input, REAL*4)

C*********************************************************************

      REAL*4 R4,RARY(IDIM*JDIM*KDIM*N)

      DO 1 K=1,IDIM*JDIM*KDIM*N
    1 RARY(K)=R4

      END


C*********************************************************************
      SUBROUTINE SETARYR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,RARY,R8)
C*********************************************************************

C  Set ALL elements of an REAL*8 grid-element array to a constant.
C  This is a work routine.

C  RARY(I,J,K) = Array to be initialized (output, REAL*8)

C  R8 = Value to which the array elements are to be set
C       (input, REAL*8)

C*********************************************************************

      REAL*8 R8,RARY(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 RARY(I,J,K)=R8

      END
C*********************************************************************
      SUBROUTINE SETARYR8N (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                      KL2,KEYOUT,NBLK,RARY,R8,N)
C*********************************************************************

C  Set ALL elements of a multidimensional REAL*8 grid-element array
C  to a constant.
C  This is a work routine.

C  RARY(I,J,K,L) = Array to be initialized (output, REAL*8)

C  R8 = Value to which the array elements are to be set
C       (input, REAL*8)

C*********************************************************************

      REAL*8 R8,RARY(IDIM*JDIM*KDIM*N)

      DO 1 K=1,IDIM*JDIM*KDIM*N
    1 RARY(K)=R8

      END

C*********************************************************************
      SUBROUTINE SETMPFAARYR8N (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &                      KL1,KL2,KEYOUT,NBLK,RARY,R8,N)
C*********************************************************************

C  Set ALL elements of a multidimensional REAL*8 grid-element array
C  to a constant.
C  This is a work routine.

C  RARY(I,J,K,L) = Array to be initialized (output, REAL*8)

C  R8 = Value to which the array elements are to be set
C       (input, REAL*8)

C*********************************************************************

      REAL*8 R8,RARY((IDIM+1)*(JDIM+1)*(KDIM+1)*N)

      DO 1 K=1,(IDIM+1)*(JDIM+1)*(KDIM+1)*N
    1 RARY(K)=R8

      END

C*********************************************************************
      SUBROUTINE SCLARYR8N (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                      KL2,KEYOUT,NBLK,RARY,R8,N)
C*********************************************************************

C  Scale ALL elements of a multidimensional REAL*8 grid-element array
C  by a constant.
C  This is a work routine.

C  RARY(I,J,K,L) = Array to be initialized (output, REAL*8)

C  R8 = Value by which the array elements are to be scaled
C       (input, REAL*8)

C*********************************************************************

      REAL*8 R8,RARY(IDIM*JDIM*KDIM*N)

      DO 1 K=1,IDIM*JDIM*KDIM*N
    1 RARY(K)=RARY(K)*R8

      END
C*********************************************************************
      SUBROUTINE CPYARYI4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,IARY1,IARY2)
C*********************************************************************

C  Copies ALL elements of one INTEGER grid-element array to another.
C  This is a work routine.

C  IARY1(I,J,K) = Source array (input, INTEGER)

C  IARY2(I,J,K) = Target array (output, INTEGER)

C*********************************************************************

      INTEGER IARY1(IDIM,JDIM,KDIM),IARY2(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 IARY2(I,J,K)=IARY1(I,J,K)

      END
C*********************************************************************
      SUBROUTINE CPYARYR4 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,RARY1,RARY2)
C*********************************************************************

C  Copies ALL elements of one REAL*4 grid-element array to another.
C  This is a work routine.

C  RARY1(I,J,K) = Source array (input, REAL*4)

C  RARY2(I,J,K) = Target array (output, REAL*4)

C*********************************************************************

      REAL*4 RARY1(IDIM,JDIM,KDIM),RARY2(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 RARY2(I,J,K)=RARY1(I,J,K)

      END
C*********************************************************************
      SUBROUTINE CPYARYR8 (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                   KL2,KEYOUT,NBLK,RARY1,RARY2)
C*********************************************************************

C  Copies ALL elements of one REAL*8 grid-element array to another.
C  This is a work routine.

C  RARY1(I,J,K) = Source array (input, REAL*8)

C  RARY2(I,J,K) = Target array (output, REAL*8)

C*********************************************************************

      REAL*8 RARY1(IDIM,JDIM,KDIM),RARY2(IDIM,JDIM,KDIM)

      DO 1 K=1,KDIM
      DO 1 J=1,JDIM
      DO 1 I=1,IDIM
    1 RARY2(I,J,K)=RARY1(I,J,K)

      END
C*********************************************************************
      SUBROUTINE CPYARYR8N(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                     KL2,KEYOUT,NBLK,RARY1,RARY2,N)
C*********************************************************************

C  Copies ALL elements of N REAL*8 grid-element arrays to another.
C  This is a work routine.

C  RARY1(I) = Source array (input, REAL*8)

C  RARY2(I) = Target array (output, REAL*8)

C*********************************************************************

      REAL*8 RARY1(IDIM*JDIM*KDIM*N),RARY2(IDIM*JDIM*KDIM*N)


      DO I = 1,IDIM*JDIM*KDIM*N
         RARY2(I) = RARY1(I)
      END DO

      END
C*********************************************************************
      SUBROUTINE COUNTE (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     &                   KEYOUT,NBLK)
C*********************************************************************

C  Counts grid elements for each physical model

C*********************************************************************
      INCLUDE 'control.h'
      INCLUDE 'layout.h'

      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)

      N=0
      DO 1 K=KL1,KL2
      JL1=JL1V(K)
      JL2=JL2V(K)
      DO 1 J=JL1,JL2
      DO 1 I=IL1,IL2
      IF (KEYOUT(I,J,K).EQ.1) N=N+1
    1 CONTINUE

      M=MODBLK(NBLK)
      NEMOD(M)=NEMOD(M)+N

! bag8 - Do not count elements twice unless driving model
!        is not the flow model.
      IF (FMODBLK(NBLK).NE.MODBLK(NBLK)) THEN
         M=FMODBLK(NBLK)
         NEMOD(M)=NEMOD(M)+N
      ENDIF

      END

C*********************************************************************
      SUBROUTINE DEPTH1_TOP (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &   KL1,KL2,KEYOUT,NBLK,DEPTH)
C*********************************************************************

C  Calculate block center depth array for orthogonal grid non-plane top
C  This is a work routine.
C  DEPTH = 0. at x=0., y=0., z=0.

C  DEPTH(I,J,K) = Block center depth in feet for fault block NBLK
C                 (output, REAL*8)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,IL1,IL2,KL1,KL2,NBLK,LDIM
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  DEPTH(IDIM,JDIM,KDIM)
      INTEGER ID,JD,KD,IG,JG,KG,I,J,K,IOFF,JOFF,KOFF,MERR
      REAL*8  DX1,DX2

      CALL BLKDIM(NBLK,ID,JD,KD,MERR)
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)

! First shift top surface down by 1/2 DX
      DO K=1,KDIM
      KG=K+KOFF
      IF (KG.GT.0.AND.KG.LE.KD) THEN
         DO J=1,JDIM
         JG=J+JOFF
         IF (JG.GT.0.AND.JG.LE.JD) THEN
            I=IL1
            IG=I+IOFF
            DX1=XREC(IG+1,NBLK)-XREC(IG,NBLK)
            DEPTH(I,J,K)=DEPTH(I,J,K) + 0.5*DX1
         ENDIF
         ENDDO
      ENDIF
      ENDDO

! Fill remaining layers
      DO K=1,KDIM
      KG=K+KOFF
      IF (KG.GT.0.AND.KG.LE.KD) THEN
         DO J=1,JDIM
         JG=J+JOFF
         IF (JG.GT.0.AND.JG.LE.JD) THEN
            DO I=IL1+1,IL2
            IG=I+IOFF
            IF (IG.GT.1.AND.IG.LE.ID) THEN
               DX1=XREC(IG,NBLK)-XREC(IG-1,NBLK)
               DX2=XREC(IG+1,NBLK)-XREC(IG,NBLK)
               DEPTH(I,J,K)=DEPTH(I-1,J,K) + 0.5*(DX1+DX2)
            ENDIF
            ENDDO
         ENDIF
         ENDDO
      ENDIF
      ENDDO

      END

C**********************************************************************
      SUBROUTINE FILL_N_MYPRC(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &   KL1,KL2,KEYOUT,NBLK,NBLK_ARR)
C**********************************************************************
C bag8 - for visualizing processor assignment
C**********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),JL2V(KDIM),
     &        KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK
      REAL*8  NBLK_ARR(IDIM,JDIM,KDIM)
      INTEGER I,J,K

      DO K = KL1,KL2
        DO J = JL1V(K),JL2V(K)
          DO I = IL1,IL2
            IF (KEYOUT(I,J,K).EQ.1) THEN
              NBLK_ARR(I,J,K) = MYPRC
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END

! --- SAUMIK,BGANIS

C*********************************************************************
      SUBROUTINE APPENDMAP(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &                     KL2,KEYOUT,NBLK,XC,YC,ZC)
C*********************************************************************
      IMPLICIT NONE
      include 'control.h'
      include 'layout.h'
      include 'blkary.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,
     &        JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM),NBLK
      REAL*8  XC(IDIM+1,JDIM+1,KDIM+1),YC(IDIM+1,JDIM+1,KDIM+1),
     &        ZC(IDIM+1,JDIM+1,KDIM+1)
      INTEGER J,K,L,JJ,KK,IOFF,JOFF,KOFF,IERR
      REAL*8  Y1,Y2,Z1,Z2,YMID,ZMID

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,IERR)

! Loop over aerial grid elements on fault block NPAY
      L=0
      DO KK = 1,NZDIM(NPAY)
      DO JJ = 1,NYDIM(NPAY)
         L=L+1
         YMID = MIDPOINTS(NPAY,1,JJ,KK)
         ZMID = MIDPOINTS(NPAY,2,JJ,KK)

! Loop over aerial grid elements on primary fault block
         DO K = KL1,KL2
         DO J = JL1V(K),JL2V(K)
            IF (KEYOUT(IL1,J,K).NE.1) CYCLE
            IF (KNDGRD.EQ.3) THEN
              Y1 = YC(IL1,J,K)
              Y2 = YC(IL1+1,J+1,K+1)
              Z1 = ZC(IL1,J,K)
              Z2 = ZC(IL1+1,J+1,K+1)
            ELSE
              Y1 = YREC(JOFF+J,NBLK)
              Y2 = YREC(JOFF+J+1,NBLK)
              Z1 = ZREC(KOFF+K,NBLK)
              Z2 = ZREC(KOFF+K+1,NBLK)
            ENDIF

! saumik - changed from .LE. to .LT.
            IF(((YMID.GE.Y1).AND.(YMID.LT.Y2)).AND.
     &         ((ZMID.GE.Z1).AND.(ZMID.LT.Z2))) THEN

               PRCMAP(START1+L) = MYPRC
               FLAGBLK(NPAY)=.TRUE.
               GO TO 1
            ENDIF
         ENDDO
         ENDDO
   1     CONTINUE
      ENDDO
      ENDDO

      END

c======================================================================
      real*8 function Tri_AREA(p1,p2,p3)
c======================================================================
      implicit none
C
c     area of triangle  1/2 |a*b| = 1/2|a||b|sin(theta),
c     a and b are edge vectors
c     a = (a1,a2,a3)  b = (b1,b2,b3)
c     |a*b|^2 = |a|^2|b|^2 - (a dot b)^2
c
c
c          p3
c         /  \
c        /    \
c       p1 --- p2
C
C
      real*8 p1(3),p2(3),p3(3)
C
      integer dim
      real*8  a(3),b(3)
C
      do 100 dim = 1,3
         a(dim) = p3(dim) - p1(dim)
         b(dim) = p2(dim) - p1(dim)
 100  continue

c  triagle p3-p1-p2
      Tri_AREA = 0.5d0*DSQRT((a(1)**2+a(2)**2+a(3)**2)*
     &     (b(1)**2+b(2)**2+b(3)**2)-
     &     (a(1)*b(1)+a(2)*b(2)+a(3)*b(3))**2)

      return
      end

C**********************************************************************
C                   END OF FILE
C**********************************************************************

