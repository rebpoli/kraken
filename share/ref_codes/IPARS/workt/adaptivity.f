!*********************************************************************
! adaptivity.df: subroutines for grid adaptivity
!*********************************************************************
! Created by Ben Ganis
! 2/28/17

C*********************************************************************
      MODULE adaptmod
C*********************************************************************
      IMPLICIT NONE
      SAVE

      INTEGER :: NLEVA       ! Number of refinement levels
      INTEGER :: NDIMA       ! Number of space dimensions
      INTEGER :: NSTRIDE     ! Number of 1d cells between levels in tree
      INTEGER :: NPHA        ! Number of phases in model
      INTEGER :: IDIM_F,JDIM_F,KDIM_F,LDIM_F  ! Local fine scale dims
      INTEGER, ALLOCATABLE :: FINELEV(:,:,:)
      REAL*8, ALLOCATABLE :: FINEARR(:,:,:,:)
      REAL*8, ALLOCATABLE :: INDICATOR(:,:,:)
      INTEGER, ALLOCATABLE :: NEWLEV(:,:,:)

! A priori adaptivity
      INTEGER, PARAMETER :: MXADAPT = 5  ! Max number of adapt vars
      INTEGER :: ADAPT_VARS  ! Number of variables on which to adapt
      CHARACTER*20 :: ADAPT_VAR(MXADAPT)      ! GEA names to adapt
      REAL*8 :: ADAPT_TOL(MXADAPT)            ! Tolerances to adapt
      INTEGER :: IVAR     ! Counter for adaptivity variables
      INTEGER :: ILEV     ! Counter for adaptivity levels

! A posteriori adaptivity
      REAL*8 :: REFINE_TOL
      REAL*8 :: COARSE_TOL

      CONTAINS

!---------------------------------------------------------------------
      SUBROUTINE UPDATE_FINEARR_I4(ARR)
!---------------------------------------------------------------------
      USE scrat1mod
      IMPLICIT NONE
      INTEGER, ALLOCATABLE :: ARR(:,:,:)
      INTEGER KTMP
      KTMP=1
      CALL I4UPDATE2(IDIM_F,JDIM_F,KDIM_F,1,1,1,NLEVA,KTMP,ARR,A)
      END SUBROUTINE UPDATE_FINEARR_I4

      END MODULE adaptmod

C*********************************************************************
      SUBROUTINE ADAPT_APRIORI()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER NERR
      EXTERNAL NEWLEVTOKEYOUT,PRINT_KEYOUT

      CALL EV_MARK_APRIORI()
      IF ((ADAPT_MODE.EQ.2).OR.(ADAPT_MODE.EQ.3)) THEN
        CALL ADAPTIVITY_SAVE_LEVEL(NERR)
        CALL CALLWORK(NEWLEVTOKEYOUT,[1,N_MARK])
!        CALL CALLWORK(PRINT_KEYOUT,0)
        CALL ADAPT_COUNT_ELE()
        IF (MYPRC.EQ.0) THEN
           WRITE(*,*)'In ADAPT_APRIORI, NEMOD=',NEMOD(13)
C           WRITE(*,*)'In ADAPT_APRIORI, NEMOD=',NEMOD(3)
        ENDIF
        CALL UPDATE(0,2)
        XAPRIORI=1
        CALL ADAPT_REINIT(NERR)
        CALL IPARS_UPDATE_PERM()
      ENDIF
      IF (ADAPT_MODE.EQ.3) THEN
C        CALL XFLASH_REINIT()
        IF (MYPRC.EQ.0) WRITE(*,*)'Calling STEP in ADAPT_APRIORI'
        CALL STEP(NERR)
      ENDIF
      END SUBROUTINE ADAPT_APRIORI

C*********************************************************************
      SUBROUTINE RESTORE_APRIORI()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER NERR
      EXTERNAL FINEKEYOUT,COUNTE

      IF ((ADAPT_MODE.EQ.2).OR.(ADAPT_MODE.EQ.3)) THEN
        CALL ADAPTIVITY_SAVE_LEVEL(NERR)
        CALL CALLWORK(FINEKEYOUT,0)
        CALL ADAPT_COUNT_ELE()
        IF (MYPRC.EQ.0) THEN
           WRITE(*,*)'In RESTORE_APRIORI, NEMOD=',NEMOD(13)
C           WRITE(*,*)'In RESTORE_APRIORI, NEMOD=',NEMOD(3)
        ENDIF
        CALL UPDATE(0,2)
        XAPRIORI=2
        CALL ADAPT_REINIT(NERR)
        CALL IPARS_UPDATE_PERM()
C        CALL XFLASH_REINIT()
      ENDIF
      END SUBROUTINE RESTORE_APRIORI

C*********************************************************************
      SUBROUTINE ADAPT_APOSTERIORI()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER N,NERR
      EXTERNAL NEWLEVTOKEYOUT

      CALL EV_MARK_APOSTERIORI()
      IF ((ADAPT_MODE.EQ.5).OR.(ADAPT_MODE.EQ.6)) THEN
        CALL ADAPTIVITY_SAVE_LEVEL(NERR)
        CALL CALLWORK(NEWLEVTOKEYOUT,[1,N_MARK])
!        CALL CALLWORK(PRINT_KEYOUT,0)
        CALL ADAPT_COUNT_ELE()
        IF (MYPRC.EQ.0) THEN
           N=NEMOD(13)
C           N=NEMOD(3)
           WRITE(*,*)'In ADAPT_APOSTERIORI, NEMOD=',N
        ENDIF
        CALL UPDATE(0,2)
        CALL ADAPT_REINIT(NERR)
        CALL IPARS_UPDATE_PERM()
C        CALL XFLASH_REINIT()
        IF (ADAPT_MODE.EQ.6) THEN
          CALL KILL_IPARS('ADAPT_MODE=6 not implemented yet')
        ENDIF
      ENDIF

      END SUBROUTINE ADAPT_APOSTERIORI

C*********************************************************************
      SUBROUTINE RESTORE_APOSTERIORI()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      END

C*********************************************************************
      SUBROUTINE ADAPT_COUNT_ELE()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      EXTERNAL COUNTE
      NEMOD(:)=0
      CALL CALLWORK(COUNTE,0)
      CALL SUMITI(19,NEMOD)
      CALL SPREAD(19,NEMOD)
      END SUBROUTINE ADAPT_COUNT_ELE

C*********************************************************************
      SUBROUTINE ADAPT_REINIT(NERR)
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER NERR

      IF (MYPRC.EQ.0) WRITE(*,*)'In ADAPT_REINIT'
      NERR = 0
      CALL GET_GELEI_LSIZE(NERR)
      IF (ADAPT_MODE.EQ.2) THEN
      CALL TUPSCALE2(NERR)
C        CALL XUPSCALE2(NERR)
      ELSEIF (ADAPT_MODE.EQ.3) THEN
      CALL TUPSCALE3(NERR)
C        CALL XUPSCALE3(NERR)
      ELSEIF ((ADAPT_MODE.EQ.4).OR.(ADAPT_MODE.EQ.5)
     &        .OR.(DEALII)) THEN
      CALL TUPSCALE1(NERR)
C        CALL XUPSCALE1(NERR)
      ELSE
        CALL KILL_IPARS('Unknown ADAPT_MODE in ADAPT_REINIT')
      ENDIF
      CALL SBLKIN(NERR)
      CALL DUALINIT(NERR)
      CALL TIMON(12)
C      MODACT=3
C      IF (MODELON(3).OR.FLOWMODEL.EQ.MODACT)
C     &             CALL XBLKIN(NERR)
C      MODACT=16
C      IF (MODELON(16).OR.FLOWMODEL.EQ.MODACT)
C     &             CALL XBLKIN(NERR)
      MODACT=0
      CALL TIMOFF(12)
      IF (NERR.GT.0) CALL KILL_IPARS('Error in DUALINIT')
      CALL HYPRE_EVFEM_INIT()

      END SUBROUTINE ADAPT_REINIT

C*********************************************************************
      SUBROUTINE IPARS_UPDATE_PERM()
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'blkary.h'
      INCLUDE 'layout.h'
      INTEGER NERR
      EXTERNAL TRANC1
C      EXTERNAL TRANC2

      IF (KNDGRD.EQ.1) THEN
        CALL UPDATE(N_XPERM,1)
        CALL UPDATE(N_YPERM,1)
        CALL UPDATE(N_ZPERM,1)
        CALL CALLWORK(TRANC1,[6,N_TCOFX,N_TCOFY,N_TCOFZ,N_XPERM,
     &    N_YPERM,N_ZPERM])
      ELSEIF (KNDGRD.EQ.3) THEN
        CALL UPDATE(N_XPERM,2)
        CALL UPDATE(N_YPERM,2)
        CALL UPDATE(N_ZPERM,2)
C        CALL CALLWORK(TRANC2,[9,N_TCOFX,N_TCOFY,N_TCOFZ,N_XPERM,
C     &    N_YPERM,N_ZPERM,N_XC,N_YC,N_ZC])
      ENDIF
      CALL IFTRAN()
      CALL SDPWELL1(NERR)

      END SUBROUTINE IPARS_UPDATE_PERM

C*********************************************************************
      SUBROUTINE EV_MARK_APRIORI()
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER  NTYPGA,ND4GA,KARR,KERR
      EXTERNAL FILL_MARK,FILL_REFINE

      IF (MYPRC.EQ.0) WRITE(*,*)'In EV_MARK_APRIORI'

      NEWLEV(:,:,:)=0
      DO IVAR=1,ADAPT_VARS
        CALL ARYDAT(ADAPT_VAR(IVAR),NTYPGA,ND4GA,KARR,KERR)
        IF ((NTYPGA.NE.2).OR.(ND4GA.NE.1)) THEN
          CALL KILL_IPARS('Type of ADAPT_VAR not supported')
        ENDIF
        CALL UPDATE(KARR,1)
        CALL CALLWORK(FILL_MARK,[2,KARR,N_MARK])
        CALL CALLWORK(FILL_REFINE,[1,N_MARK])
      ENDDO
      CALL EV_GRADE_NEWLEV()

      END SUBROUTINE EV_MARK_APRIORI

C*********************************************************************
      SUBROUTINE EV_MARK_APOSTERIORI()
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'blkary.h'
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'tarydat.h'
C      INCLUDE 'xarydat.h'
      INTEGER  NTYPGA,ND4GA,KARR,KERR,NRHS,NVEL
      EXTERNAL FILL_MARK_APOSTERIORI_RESIDUAL,
     &  FILL_MARK_APOSTERIORI_CURVATURE,FILL_REFINE,FILL_COARSEN,
     &  KEYOUT_NEWLEV0
      LOGICAL :: ONCEONLY = .TRUE.

      IF (MYPRC.EQ.0) WRITE(*,*)'In EV_MARK_APOSTERIORI'

! Reset new level information
      IF (ADAPT_RESET) THEN
        NEWLEV(:,:,:)=1
      ELSEIF (ONCEONLY) THEN
        CALL CALLWORK(KEYOUT_NEWLEV0,0)
        ONCEONLY = .FALSE.
      ENDIF

! Compute weighted residual
      NRHS=N_RESIDV(13)
      NVEL=N_VEL
C      NRHS=N_RESIDV(3)
C      NVEL=N_XVEL
      IF (ADAPT_INDICATOR.EQ.1) THEN
        IF (NRHS.EQ.0) CALL KILL_IPARS('EV_MARK_APOSTERIORI: NRHS=0')
        CALL CALLWORK(FILL_MARK_APOSTERIORI_RESIDUAL,
     &    [2,NRHS,N_MARK])
      ELSEIF (ADAPT_INDICATOR.EQ.2) THEN
        IF (NVEL.EQ.0) CALL KILL_IPARS('EV_MARK_APOSTERIORI: NVEL=0')
        CALL CALLWORK(FILL_MARK_APOSTERIORI_CURVATURE,
     &    [2,NVEL,N_MARK])
      ELSE
        CALL KILL_IPARS('EV_MARK_APOSTERIORI: unknown indicator')
      ENDIF

! Set level where residual exceeds tolerance and grade mesh
      IF (ADAPT_MODE.NE.4) THEN
        CALL CALLWORK(FILL_REFINE,[1,N_MARK])
        CALL CALLWORK(FILL_COARSEN,[1,N_MARK])
        CALL EV_GRADE_NEWLEV()
      ENDIF

      END SUBROUTINE EV_MARK_APOSTERIORI

C*********************************************************************
      SUBROUTINE EV_GRADE_NEWLEV()
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER  NTYPGA,ND4GA,KARR,KERR
      EXTERNAL FILL_NEWLEV1,FILL_NEWLEV2,FILL_NEWLEV3,FILL_NEWLEV4

      CALL CALLWORK(FILL_NEWLEV1,[1,N_MARK])
      CALL UPDATE_FINEARR_I4(NEWLEV)
      DO ILEV=NLEVA-1,2,-1
        CALL CALLWORK(FILL_NEWLEV2,[1,N_MARK])
        CALL UPDATE_FINEARR_I4(NEWLEV)
      ENDDO
      CALL CALLWORK(FILL_NEWLEV3,[1,N_MARK])

      END SUBROUTINE EV_GRADE_NEWLEV

C*********************************************************************
      SUBROUTINE FILL_MARK(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK,VAR,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  VAR(IDIM,JDIM,KDIM),MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2
      REAL*8  DELTA

      MARK = 0.D0
      IF (NBLK.NE.NLEVA) RETURN

      DO K=KL1,KL2
      DO J=JL1V(K),JL2V(K)
      DO I=IL1,IL2
      IF (KEYOUT(I,J,K).EQ.1) THEN
        IF (ABS(KEYOUT(I-1,J,K)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I-1,J,K))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
        IF (ABS(KEYOUT(I+1,J,K)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I+1,J,K))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
        IF (ABS(KEYOUT(I,J-1,K)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I,J-1,K))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
        IF (ABS(KEYOUT(I,J+1,K)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I,J+1,K))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
        IF (ABS(KEYOUT(I,J,K-1)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I,J,K-1))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
        IF (ABS(KEYOUT(I,J,K+1)).EQ.1) THEN
          DELTA = ABS(VAR(I,J,K)-VAR(I,J,K+1))
          MARK(I,J,K)=MAX(MARK(I,J,K),DELTA)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_MARK

C*********************************************************************
      SUBROUTINE FILL_MARK_APOSTERIORI_RESIDUAL(IDIM,JDIM,KDIM,
     &   LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,RESID,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'tbaldat.h'
C      INCLUDE 'xbaldat.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  RESID(IDIM,JDIM,KDIM),MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2,IOFF,JOFF,KOFF,MERR,STRIDE,
     &  IFINE,JFINE,KFINE,ILAST,JLAST,KLAST
      REAL*8  DX,DY,DZ,H,R2

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      R2 = 1.D0
      R2 = RESID2
C      R2 = RESID2(1)
      STRIDE=NSTRIDE**(NLEVA-NBLK)
      MARK = 0.D0
      DO K=KL1,KL2
      DZ=DZREC(K+KOFF,NBLK)
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J=JL1V(K),JL2V(K)
      DY=DYREC(J+JOFF,NBLK)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I=IL1,IL2
      DX=DXREC(I+IOFF,NBLK)
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        H = SQRT(DX**2+DY**2+DZ**2)
        MARK(I,J,K)=SQRT(H**2*RESID(I,J,K)**2)/R2
        INDICATOR(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=MARK(I,J,K)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_MARK_APOSTERIORI_RESIDUAL

C*********************************************************************
      SUBROUTINE FILL_MARK_APOSTERIORI_CURVATURE(IDIM,JDIM,KDIM,
     &   LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,VEL,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  VEL(IDIM,JDIM,KDIM,3,NPHA),MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2,IOFF,JOFF,KOFF,MERR,STRIDE,
     &  IFINE,JFINE,KFINE,ILAST,JLAST,KLAST,IPH
      REAL*8  DX,DY,DZ,DUX,DUY,DUZ

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      STRIDE=NSTRIDE**(NLEVA-NBLK)
      MARK = 0.D0
      DO IPH=1,NPHA
      DO K=KL1,KL2
      DZ=DZREC(K+KOFF,NBLK)
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J=JL1V(K),JL2V(K)
      DY=DYREC(J+JOFF,NBLK)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I=IL1,IL2
      DX=DXREC(I+IOFF,NBLK)
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        DUX=ABS(VEL(I-1,J,K,1,IPH)-VEL(I,J,K,1,IPH))/DX
        DUY=ABS(VEL(I,J-1,K,2,IPH)-VEL(I,J,K,2,IPH))/DY
        DUZ=ABS(VEL(I,J,K-1,3,IPH)-VEL(I,J,K,3,IPH))/DZ
        MARK(I,J,K)=MAX(MARK(I,J,K),DUX,DUY,DUZ)
        INDICATOR(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=MARK(I,J,K)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_MARK_APOSTERIORI_CURVATURE

C*********************************************************************
      SUBROUTINE FILL_REFINE(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,STRIDE,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST

      IF ((ADAPT_MODE.GE.1).AND.(ADAPT_MODE.LE.3)) THEN
        IF (NBLK.NE.NLEVA) RETURN
      ENDIF

!----------------------------------------------------------------------
! Set finest level where threshold is exceeded (a priori)
! or refine by one level (a posteriori)
!----------------------------------------------------------------------
      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K=KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J=JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I=IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        IF ((ADAPT_MODE.GE.1).AND.(ADAPT_MODE.LE.3)) THEN
          IF (MARK(I,J,K).GE.ADAPT_TOL(IVAR)) THEN
            NEWLEV(I,J,K)=NLEVA
          ENDIF
        ELSEIF ((ADAPT_MODE.GE.4).AND.(ADAPT_MODE.LE.6)) THEN
          IF (MARK(I,J,K).GE.REFINE_TOL) THEN
!            WRITE(*,*)'SETTING NEWLEV, NBLK,I,J,K,MARK,TOL=',
!     &        NBLK,I,J,K,MARK(I,J,K),REFINE_TOL
            IF (ADAPT_REFINE_MOVE1) THEN
              NEWLEV(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=
     &          MIN(NBLK+1,NLEVA)
            ELSE
              NEWLEV(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=NLEVA
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_REFINE

C*********************************************************************
      SUBROUTINE FILL_COARSEN(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &  JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,STRIDE,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST,I2,J2,K2

      IF (ADAPT_MODE.LE.4) RETURN

!----------------------------------------------------------------------
! Coarsen by one level if indicator is below coarsening threshold
!----------------------------------------------------------------------
      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K=KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J=JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I=IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
! Check for an appropriate candidate element to coarsen
      IF (KEYOUT(I,J,K).EQ.1) GOTO 10
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        IF ((ADAPT_REFINE_MOVE1).AND.(FINELEV(I2,J2,K2).GT.NBLK+1))
     &    GOTO 10
        IF (FINELEV(I2,J2,K2).LE.NBLK) GOTO 10
      ENDDO
      ENDDO
      ENDDO
! Coarsen if entire indicator region is below threshold
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        IF (INDICATOR(I2,J2,K2).GT.COARSE_TOL) GOTO 10
      ENDDO
      ENDDO
      ENDDO
      NEWLEV(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=NBLK
 10   CONTINUE
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_COARSEN

C*********************************************************************
      SUBROUTINE FILL_NEWLEV1(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INCLUDE 'wells.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2,STRIDE,
     &        IMAX,JMAX,KMAX,IFINE,JFINE,KFINE,LEV,ILAST,JLAST,KLAST,
     &        IPREV,JPREV,KPREV,INEXT,JNEXT,KNEXT,I2,J2,K2,
     &        DIR,MAXDIR,IFINE2,JFINE2,KFINE2,ILAST2,JLAST2,KLAST2,
     &        IW,L,IOFF,JOFF,KOFF,MERR
      LOGICAL FILL

      IF (NBLK.NE.NLEVA) RETURN

!----------------------------------------------------------------------
! Also set finest level where wells were present
!----------------------------------------------------------------------

      MODACT=13
C      MODACT=3
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      DO IW=1,NUMWEL
      IF (MODWEL(IW).NE.MODACT
C     & .AND.MODACT.NE.FLOWMODEL
     &   ) CYCLE
      DO L=1,NUMELE(IW)
      IF (LOCWEL(6,L,IW).EQ.MYPRC) THEN
        IF (LOCWEL(1,L,IW).NE.NLEVA) THEN
          CALL KILL_IPARS('Wells need to start on finest grid level')
        ENDIF
        I=LOCWEL(3,L,IW)-IOFF
        J=LOCWEL(4,L,IW)-JOFF
        K=LOCWEL(5,L,IW)-KOFF
        NEWLEV(I,J,K)=NLEVA
      ENDIF
      ENDDO
      ENDDO
      MODACT=0

!----------------------------------------------------------------------
! Next make finest level consistent
!----------------------------------------------------------------------
      STRIDE=NSTRIDE
      KMAX=(KDIM-2*KLAY)/STRIDE+KLAY
      JMAX=(JDIM-2*JLAY)/STRIDE+JLAY
      IF (NDIMA.EQ.3) THEN
        IMAX=(IDIM-2*ILAY)/STRIDE+ILAY
      ELSE
        IMAX=IDIM-ILAY
      ENDIF
      DO K = KLAY+1,KMAX
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JLAY+1,JMAX
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = ILAY+1,IMAX
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF

      FILL=.FALSE.
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        IF (NEWLEV(I2,J2,K2).EQ.NLEVA) THEN
          FILL=.TRUE.
          GOTO 10
        ENDIF
      ENDDO
      ENDDO
      ENDDO

 10   CONTINUE
      IF (FILL) THEN
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
          NEWLEV(I2,J2,K2)=MAX(NEWLEV(I2,J2,K2),NLEVA)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_NEWLEV1

C*********************************************************************
      SUBROUTINE FILL_NEWLEV2(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2,STRIDE,
     &        IMAX,JMAX,KMAX,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST,
     &        IPREV,JPREV,KPREV,INEXT,JNEXT,KNEXT,I2,J2,K2,
     &        DIR,MAXDIR,IFINE2,JFINE2,KFINE2,ILAST2,JLAST2,KLAST2
      LOGICAL FILL

      IF (NBLK.NE.NLEVA) RETURN

!---------------------------------------------------------------------
! Loop over coarser levels, grading the mesh
!----------------------------------------------------------------------
      STRIDE=NSTRIDE**(NLEVA-ILEV+1)
      KMAX=(KDIM-2*KLAY)/STRIDE+KLAY
      JMAX=(JDIM-2*JLAY)/STRIDE+JLAY
      IF (NDIMA.EQ.3) THEN
        IMAX=(IDIM-2*ILAY)/STRIDE+ILAY
      ELSE
        IMAX=IDIM-ILAY
      ENDIF

      !----------------------------------------------------------------
      ! Check to see if there are any refinements > ILEV in
      ! neighboring blocks.  If so, fill my ILEV block with ILEV.
      !----------------------------------------------------------------
      DO K = KLAY+1,KMAX
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JLAY+1,JMAX
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = ILAY+1,IMAX
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF

      IF (NDIMA.EQ.3) THEN
        MAXDIR=6
      ELSE
        MAXDIR=4
      ENDIF

      DO DIR=1,MAXDIR

      KFINE2=KFINE
      KLAST2=KLAST
      JFINE2=JFINE
      JLAST2=JLAST
      IFINE2=IFINE
      ILAST2=ILAST
      IF (DIR.EQ.1) THEN
        JFINE2=JFINE-1
        JLAST2=JFINE2
      ELSEIF (DIR.EQ.2) THEN
        JFINE2=JLAST+1
        JLAST2=JFINE2
      ELSEIF (DIR.EQ.3) THEN
        KFINE2=KFINE-1
        KLAST2=KFINE2
      ELSEIF (DIR.EQ.4) THEN
        KFINE2=KLAST+1
        KLAST2=KFINE2
      ELSEIF (DIR.EQ.5) THEN
        IFINE2=IFINE-1
        ILAST2=IFINE2
      ELSE
        IFINE2=ILAST+1
        ILAST2=IFINE2
      ENDIF

      FILL=.FALSE.
      DO K2=KFINE2,KLAST2
      DO J2=JFINE2,JLAST2
      DO I2=IFINE2,ILAST2
        IF (NEWLEV(I2,J2,K2).GT.ILEV) THEN
          FILL=.TRUE.
          GOTO 40
        ENDIF
      ENDDO
      ENDDO
      ENDDO

 40   CONTINUE
      IF (FILL) THEN
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        NEWLEV(I2,J2,K2)=MAX(NEWLEV(I2,J2,K2),ILEV)
      ENDDO
      ENDDO
      ENDDO
      GOTO 50
      ENDIF

      ENDDO ! End loop over DIR

 50   CONTINUE
      ENDDO ! End loop over I
      ENDDO ! End loop over J
      ENDDO ! End loop over K

      !----------------------------------------------------------------
      ! Check to see if there are any refinements > ILEV in my ILEV
      ! block. If so, fill the remainder of the ILEV block with ILEV.
      !----------------------------------------------------------------
      DO K = KLAY+1,KMAX
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JLAY+1,JMAX
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = ILAY+1,IMAX
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF

      FILL=.FALSE.
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        IF (NEWLEV(I2,J2,K2).GT.ILEV) THEN
          FILL=.TRUE.
          GOTO 20
        ENDIF
      ENDDO
      ENDDO
      ENDDO
 20   CONTINUE
      IF (FILL) THEN
      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        NEWLEV(I2,J2,K2)=MAX(NEWLEV(I2,J2,K2),ILEV)
      ENDDO
      ENDDO
      ENDDO
      GOTO 30
      ENDIF

 30   CONTINUE
      ENDDO ! End loop over I
      ENDDO ! End loop over J
      ENDDO ! End loop over K

      END SUBROUTINE FILL_NEWLEV2

C*********************************************************************
      SUBROUTINE FILL_NEWLEV3(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8  MARK(IDIM,JDIM,KDIM)
      INTEGER I,J,K,JL1,JL2,STRIDE,
     &        IMAX,JMAX,KMAX,IFINE,JFINE,KFINE,LEV,ILAST,JLAST,KLAST,
     &        IPREV,JPREV,KPREV,INEXT,JNEXT,KNEXT,I2,J2,K2,
     &        DIR,MAXDIR,IFINE2,JFINE2,KFINE2,ILAST2,JLAST2,KLAST2
      LOGICAL FILL

      IF (NBLK.NE.NLEVA) RETURN

!----------------------------------------------------------------------
! Fill any remainder with coarsest level
!----------------------------------------------------------------------
      LEV=1
      STRIDE=NSTRIDE**(NLEVA-LEV)
      KMAX=(KDIM-2*KLAY)/STRIDE+KLAY
      JMAX=(JDIM-2*JLAY)/STRIDE+JLAY
      IF (NDIMA.EQ.3) THEN
        IMAX=(IDIM-2*ILAY)/STRIDE+ILAY
      ELSE
        IMAX=IDIM-ILAY
      ENDIF
      DO K = KLAY+1,KMAX
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JLAY+1,JMAX
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = ILAY+1,IMAX
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF

      DO K2=KFINE,KLAST
      DO J2=JFINE,JLAST
      DO I2=IFINE,ILAST
        IF (NEWLEV(I2,J2,K2).EQ.0) THEN
          NEWLEV(I2,J2,K2)=LEV
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE FILL_NEWLEV3

C*********************************************************************
      SUBROUTINE NEWLEVTOKEYOUT(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &     KL1,KL2,KEYOUT,NBLK,MARK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,IFINE,JFINE,KFINE,STRIDE
      REAL*8  MARK(IDIM,JDIM,KDIM)

      KEYOUT(:,:,:)=0

      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K = KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      DO J = JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      DO I = IL1,IL2
      IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
      IF (NEWLEV(IFINE,JFINE,KFINE).EQ.NBLK) THEN
        KEYOUT(I,J,K)=1
        IF (ADAPT_MARK_LEVEL) MARK(I,J,K)=DFLOAT(NBLK)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END

C*********************************************************************
      SUBROUTINE FINEKEYOUT(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &     KL1,KL2,KEYOUT,NBLK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,IFINE,JFINE,KFINE,STRIDE

      KEYOUT(:,:,:)=0
      IF (NBLK.NE.NLEVA) RETURN

      ! Set KEYOUT=1 for NBLK=NLEVA (Finest level)
      DO K = KL1,KL2
      DO J = JL1V(K),JL2V(K)
      DO I = IL1,IL2
        KEYOUT(I,J,K)=1
      ENDDO
      ENDDO
      ENDDO

      END

C*********************************************************************
      SUBROUTINE EV_COARSE_KEYOUT()
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER  NERR
      EXTERNAL COARSEKEYOUT,NEWLEVTOKEYOUT,PRINT_KEYOUT

      IF (MYPRC.EQ.0) WRITE(*,*)'In EV_COARSE_KEYOUT'

! First set coarse grid everywhere
      CALL CALLWORK(COARSEKEYOUT,0)

! Make sure well elements are on fine grid, and grade the mesh
      NEWLEV(:,:,:)=0
      CALL EV_GRADE_NEWLEV()
      CALL ADAPTIVITY_SAVE_LEVEL(NERR)
      CALL CALLWORK(NEWLEVTOKEYOUT,[1,N_MARK])
!      CALL CALLWORK(PRINT_KEYOUT,0)
      CALL ADAPT_COUNT_ELE()
      CALL UPDATE(0,2)
      CALL ADAPT_REINIT(NERR)
      CALL IPARS_UPDATE_PERM()
C        CALL XFLASH_REINIT()

      END SUBROUTINE EV_COARSE_KEYOUT

C*********************************************************************
      SUBROUTINE COARSEKEYOUT(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &     KL1,KL2,KEYOUT,NBLK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,IFINE,JFINE,KFINE,STRIDE

      KEYOUT(:,:,:)=0
      IF (NBLK.NE.1) RETURN

      ! Set KEYOUT=1 for NBLK=1 (Coarsest level)
      DO K = KL1,KL2
      DO J = JL1V(K),JL2V(K)
      DO I = IL1,IL2
        KEYOUT(I,J,K)=1
      ENDDO
      ENDDO
      ENDDO

      END

C*********************************************************************
      SUBROUTINE GETFINEDIM(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,
     &     KL2,KEYOUT,NBLK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
C      INCLUDE 'xmodel.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER IERR
      IF (NBLK.EQ.NUMBLK) THEN
        IDIM_F=IDIM
        JDIM_F=JDIM
        KDIM_F=KDIM
        LDIM_F=1
        NPHA=1
C        LDIM_F=MAX(NPH,NC,NAQ)
C        NPHA=NPH
      ENDIF
      END

C*********************************************************************
      SUBROUTINE ADAPTIVITY_SAVE_LEVEL(IERR)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INTEGER IERR
      EXTERNAL SAVE_LEVEL
      FINELEV(:,:,:)=0
      CALL CALLWORK(SAVE_LEVEL,0)
      END SUBROUTINE ADAPTIVITY_SAVE_LEVEL

C*********************************************************************
      SUBROUTINE KEYOUT_NEWLEV0(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,STRIDE,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST

      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K = KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        NEWLEV(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=NBLK
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE KEYOUT_NEWLEV0

C*********************************************************************
      SUBROUTINE SAVE_LEVEL(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,STRIDE,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST

      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K = KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        FINELEV(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST)=NBLK
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE SAVE_LEVEL

C*********************************************************************
      SUBROUTINE PRINT_KEYOUT(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER I,J,K,L

      WRITE(*,'(A,I1,A)')'KEYOUT',NBLK,'()='
      L=0
      DO K = KL1,KL2
      DO J = JL1V(K),JL2V(K)
      DO I = IL1,IL2
        WRITE(*,'(1X,I1)',ADVANCE='no') KEYOUT(I,J,K)
        L=L+1
        IF (L.EQ.30) THEN
          WRITE(*,*)
          L=0
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      IF (L.NE.0) WRITE(*,*)

      END SUBROUTINE PRINT_KEYOUT

C*********************************************************************
      SUBROUTINE SAVE_ARRAY(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,ARRAY,L4)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,L4
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 ARRAY(IDIM,JDIM,KDIM,L4)
      INTEGER I,J,K,L,STRIDE,IFINE,JFINE,KFINE,ILAST,JLAST,KLAST

      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO L = 1,L4
      DO K = KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (FINELEV(IFINE,JFINE,KFINE).EQ.NBLK) THEN
        FINEARR(IFINE:ILAST,JFINE:JLAST,KFINE:KLAST,L)=ARRAY(I,J,K,L)
!        WRITE(*,'(a,4i4,a,f11.4,a,3i4)')
!     &    'SAVE_ARRAY: I,J,K,NBLK=',I,J,K,NBLK,
!     &    ', FARR=',FINEARR(IFINE,JFINE,KFINE),
!     &    ', IFINE,JFINE,KFINE=',IFINE,JFINE,KFINE
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE SAVE_ARRAY

C*********************************************************************
      SUBROUTINE UPSCALE_ARRAY(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,ARRAY,L4)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,L4
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 ARRAY(IDIM,JDIM,KDIM,L4)
      INTEGER I,J,K,L,IFINE,JFINE,KFINE,STRIDE,I2,J2,K2,S2,LEV,
     &  ILAST,JLAST,KLAST,N
      REAL*8 WEIGHT,WSUM,AVG

      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO L = 1,L4
      DO K = KL1,KL2
      KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
      KLAST=(K-KLAY)*STRIDE + KLAY
      DO J = JL1V(K),JL2V(K)
      JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
      JLAST=(J-JLAY)*STRIDE + JLAY
      DO I = IL1,IL2
      IF (NDIMA.EQ.3) THEN
        IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
        ILAST=(I-ILAY)*STRIDE + ILAY
      ELSE
        IFINE=I
        ILAST=IFINE
      ENDIF
      IF (KEYOUT(I,J,K).EQ.1) THEN
        LEV=FINELEV(IFINE,JFINE,KFINE)
        IF (NBLK.EQ.LEV) THEN
          ! No upscaling necessary
        ELSEIF (NBLK.LT.LEV) THEN
          ! Upscaling case
!          WRITE(*,'(a,5i4)')'Up: IFINE,JFINE,KFINE,NBLK<LEV=',
!     &      IFINE,JFINE,KFINE,NBLK,LEV
          AVG=0.D0
          WSUM=0.D0
          WEIGHT=1.D0/((KLAST-KFINE+1)*(JLAST-JFINE+1)*(ILAST-IFINE+1))
          DO K2=KFINE,KLAST
          DO J2=JFINE,JLAST
          DO I2=IFINE,ILAST
            IF (FINELEV(I2,J2,K2).NE.0) THEN
!              WEIGHT=2.D0**(NDIMA*(NBLK-LEV))
              WSUM=WSUM+WEIGHT
              AVG=AVG+WEIGHT*FINEARR(I2,J2,K2,L)
!              WRITE(*,'(a,4i4,2e11.4)')'I2,J2,K2,LEV,FARR,W=',
!     &          I2,J2,K2,LEV,FINEARR(I2,J2,K2,L),WEIGHT
            ELSE
              WRITE(*,*)'NBLK,I,J,K,LEV=',NBLK,I,J,K,LEV
              WRITE(*,*)'WEIGHT=',WEIGHT
              WRITE(*,*)'FINELEV='
              N=JLAST-JFINE+1
              WRITE(*,'(<N>I2)')FINELEV(IFINE,JFINE:JLAST,KFINE:KLAST)
              CALL KILL_IPARS('error: lev=0 during upscaling')
            ENDIF
          ENDDO
          ENDDO
          ENDDO
          IF (ABS(WSUM-1.D0).GT.1.D-10) THEN
            WRITE(*,*)'NBLK,I,J,K,LEV=',NBLK,I,J,K,LEV
            WRITE(*,*)'WEIGHT=',WEIGHT
            WRITE(*,*)'FINELEV='
            N=JLAST-JFINE+1
            WRITE(*,'(<N>I2)')FINELEV(IFINE,JFINE:JLAST,KFINE:KLAST)
            WRITE(*,*)'AVG,WSUM=',AVG,WSUM
            CALL KILL_IPARS('error: upscaling failed, wsum not 1')
          ENDIF
          ARRAY(I,J,K,L)=AVG
        ELSE
          ! Downscaling case
!          WRITE(*,'(a,6i4)')'Down: MYPRC,IFINE,JFINE,KFINE,NBLK>LEV=',
!     &      MYPRC,IFINE,JFINE,KFINE,NBLK,LEV
          ARRAY(I,J,K,L)=FINEARR(IFINE,JFINE,KFINE,L)
        ENDIF
 10   ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE UPSCALE_ARRAY

C*********************************************************************
      SUBROUTINE PVTOPOR(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &                   KL1,KL2,KEYOUT,NBLK,PV)
C*********************************************************************
      USE adaptmod
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 PV(IDIM,JDIM,KDIM)
      REAL*4 DX,DY,DZ
      INTEGER I,J,K,MERR,IOFF,JOFF,KOFF,STRIDE,KFINE,JFINE,IFINE

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      STRIDE=NSTRIDE**(NLEVA-NBLK)
      DO K=KL1,KL2
         KFINE=(K-KLAY-1)*STRIDE + 1 + KLAY
         DZ=DZREC(K+KOFF,NBLK)
         DO J=JL1V(K),JL2V(K)
            JFINE=(J-JLAY-1)*STRIDE + 1 + JLAY
            DY=DYREC(J+JOFF,NBLK)
            DO I=IL1,IL2
               IFINE=(I-ILAY-1)*STRIDE + 1 + ILAY
               DX=DXREC(I+IOFF,NBLK)
               IF (FINELEV(IFINE,JFINE,KFINE).EQ.NBLK) THEN
                  PV(I,J,K)=PV(I,J,K)/(DX*DY*DZ)
!                  WRITE(*,'(A,4I4,F9.2)')'PVTOPOR: NBLK,I,J,K,POR=',
!     &              NBLK,I,J,K,PV(I,J,K)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      END

C*********************************************************************
      SUBROUTINE PORTOPV(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,
     &                   KL1,KL2,KEYOUT,NBLK,PV)
C*********************************************************************
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'layout.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 PV(IDIM,JDIM,KDIM)
      REAL*4 DX,DY,DZ
      INTEGER I,J,K,MERR,IOFF,JOFF,KOFF

      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,MERR)
      DO K=KL1,KL2
         DZ=DZREC(K+KOFF,NBLK)
         DO J=JL1V(K),JL2V(K)
            DY=DYREC(J+JOFF,NBLK)
            DO I=IL1,IL2
               IF (KEYOUT(I,J,K).EQ.1) THEN
                  DX=DXREC(I+IOFF,NBLK)
                  PV(I,J,K)=PV(I,J,K)*DX*DY*DZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      END

