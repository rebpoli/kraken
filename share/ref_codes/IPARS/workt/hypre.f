C  HYPRE.F - HYPRE INTERFACE
C
C  ROUTINES IN THIS MODULE:
C
C SUBROUTINE HYPREI(NERR,KINP)
C SUBROUTINE HYPRE_ALLOCATE(KERR)
C SUBROUTINE GET_GELEI_LSIZE(NERR)
C SUBROUTINE CALC_LSIZE(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
C                JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
C SUBROUTINE CALC_GELEI(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
C                JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,GELEI)
C SUBROUTINE CALC_ILOW_IHIGH(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
C                JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,GELEI)
C SUBROUTINE GETCOLSVALS_7PT(I,J,K,COLS,VALUES,NNZ,KEYOUT,
C                GELEI,COF,IDIM,JDIM,KDIM)
C SUBROUTINE GETCOLSVALS_7PT_R8(I,J,K,COLS,VALUES,NNZ,KEYOUT,
C                GELEI,COF,IDIM,JDIM,KDIM)
C SUBROUTINE GETCOLSVALS_7PT_2VAR(I,J,K,COLS,VALUES,NNZ,KEYOUT,
C                GELEI,COF,IDIM,JDIM,KDIM)
C SUBROUTINE GETCOLSVALS_27PT(I,J,K,COLS,VALUES,NNZ,KEYOUT,
C                GELEI,COF,IDIM,JDIM,KDIM)
C SUBROUTINE HYPRE_SOLVE(ITLIN,NERR)
C SUBROUTINE HYPRE_TRCHEM(ITLIN,NERR)
C SUBROUTINE HYPRE_IPARS_EVFEM(COFARR,,DUNKARR,RESIDARR,ITLIN)
C SUBROUTINE HYPRE_EVFEM_FILL(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
C    JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,COF,RESID,GELEI,ROWS,WORKSOL)
C SUBROUTINE HYPRE_EVFEM_SOLN(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
C    JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,DUNK,GELEI,ROWS,WORKSOL)

C======================================================================
      SUBROUTINE HYPREI(NERR,KINP)
C======================================================================

CGUS ROUTINE INITIALIZES HYPRE SOLVER

C NERR = Error number stepped by 1 on error (input & output, INTEGER)

C KINP = Input
C      = 1 ==> initial data
C      = 2 ==> transient data

C======================================================================
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'hypre.h'
      INTEGER NERR,KINP,NDUM
C
C HYPRE SOLVER PARAMETERS
C
      DBGEV=.FALSE.            ! bag8 - HYPRE_SOLVE debugging

      IF (KINP==1) THEN
        LSOL_TOL     = 1.D-9  ! Linear solver relative tolerance
        LSOL_ITMAX   = 100    ! Maximum linear solver iterations
        HYPRE_SOL_ID = 2      ! 0=AMG, 1=CG+AMG, 2=GMRES+AMG
        PRECOND_ID   = 2
        K_DIM        = 50     ! Max Krylov subspace dim before restart
        RELAX_TYPE   = 6
        CYCLE_TYPE   = 1      ! 1=V-cycle, 2=W-cycle
        COARSEN_TYPE = 6      ! 6=Falgout coarsening.
                              ! Others: 0,1,3,7,8,9,10,11.  See hypre manual.
        MEASURE_TYPE = 0
        SETUP_TYPE   = 1
        NUM_SWEEPS   = 1
        MAX_LEVELS   = 20
        TRUNC_FACTOR = 0.D0
        STRONG_THRES = 0.95D0

        PRINT_SOL    = 0      ! PRINT_SOL = 0,1
        PRINT_LEVEL  = 0      ! PRINT_LEVEL = 0,1,2,3
        HYPRE_EVFEM = .FALSE.
      ENDIF

      CALL GETVAL('LSOL_TOL ',LSOL_TOL,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,1) LSOL_TOL
      ENDIF
    1 FORMAT('Linear solver relative tolerance =',T45,G12.4)

      CALL GETVAL('LSOL_ITMAX ',LSOL_ITMAX,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,2) LSOL_ITMAX
      ENDIF
    2 FORMAT('Linear solver max. iterations =',T45,I5)

      CALL GETVAL('HYPRE_SOL_ID ',HYPRE_SOL_ID,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,3) HYPRE_SOL_ID
      ENDIF
    3 FORMAT('Hypres solver and precond. selections  =',T45,I5)

      CALL GETVAL('PRECOND_ID ',PRECOND_ID,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,17) PRECOND_ID
      ENDIF
   17 FORMAT('Hypre preconditioner ID  =',T45,I5)

      IF (HYPRE_SOL_ID/=2) STRONG_THRES = 0.25D0

      CALL GETVAL('RELAX_TYPE ',RELAX_TYPE,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,4) RELAX_TYPE
      ENDIF
    4 FORMAT('Relaxation type =',T45,I5)

      CALL GETVAL('CYCLE_TYPE ',CYCLE_TYPE,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,5) CYCLE_TYPE
      ENDIF
    5 FORMAT('Cycle type =',T45,I5)

      CALL GETVAL('COARSEN_TYPE ',COARSEN_TYPE,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,12) COARSEN_TYPE
      ENDIF
   12 FORMAT('COARSEN_TYPE =',T45,I5)

      CALL GETVAL('K_DIM ',K_DIM,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,13) K_DIM
      ENDIF
   13 FORMAT('K_DIM =',T45,I5)

      CALL GETVAL('NUM_SWEEPS ',NUM_SWEEPS,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,14) NUM_SWEEPS
      ENDIF
   14 FORMAT('NUM_SWEEPS =',T45,I5)

      CALL GETVAL('MAX_LEVELS ',MAX_LEVELS,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,15) MAX_LEVELS
      ENDIF
   15 FORMAT('MAX_LEVELS =',T45,I5)

      CALL GETVAL('PRINT_SOL ',PRINT_SOL,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,6) PRINT_SOL
      ENDIF
    6 FORMAT('Print hypre solution =',T45,I5)

      CALL GETVAL('PRINT_LEVEL ',PRINT_LEVEL,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,7) PRINT_LEVEL
      ENDIF
    7 FORMAT('Print level for hypre =',T45,I5)

      CALL GETVAL('MEASURE_TYPE ',MEASURE_TYPE,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,8) MEASURE_TYPE
      ENDIF
    8 FORMAT('Measure type for hypre =',T45,I5)

      CALL GETVAL('SETUP_TYPE ',SETUP_TYPE,'I4',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,9) SETUP_TYPE
      ENDIF
    9 FORMAT('Setup type for hypre =',T45,I5)

      CALL GETVAL('STRONG_THRES ',STRONG_THRES,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,10) STRONG_THRES
      ENDIF
   10 FORMAT('Strong threshold for hypre =',T45,G12.4)

      CALL GETVAL('TRUNC_FACTOR ',TRUNC_FACTOR,'R8',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,11) TRUNC_FACTOR
      ENDIF
   11 FORMAT('Truncation factor for hypre =',T45,G12.4)

      CALL GETVAL('HYPRE_EVFEM ',HYPRE_EVFEM,'FG',0,0,0,0,NDUM,NERR)
      IF (LEVELC.AND.((NDUM>0).OR.(KINP==1))) THEN
         WRITE(NFOUT,16) HYPRE_EVFEM
      ENDIF
   16 FORMAT('HYPRE_EVFEM fully coupled flag =',T48,L)

      END SUBROUTINE HYPREI

C======================================================================
      SUBROUTINE HYPRE_ALLOCATE(KERR)
C======================================================================
C
C ALLOCATE SPACE FOR HYPRE
C
      IMPLICIT NONE
      INCLUDE 'hypre.h'
      INCLUDE 'control.h'

      INTEGER KERR,KE
      LOGICAL :: DBG = .FALSE.

      KE = 0
      CALL ALCGEA('GELEINDEX ',4,0,N_GELEI,KERR)
      IF (KERR.GT.0) GO TO 13
      KE = KE + 1
      CALL ALCGEA ('HYPREIROW ',4,0,N_HYPRE_ROWS,KERR)
      IF (KERR.GT.0) GO TO 13
      KE = KE + 1
      CALL ALCGEA ('HYPREWORK ',2,0,N_HYPRE_WORK,KERR)
      IF (KERR.GT.0) GO TO 13

      IF (DBG) THEN
         WRITE(0,'(A, 3(A,I4))') 'HYPRE_ALLOCATE: ','N_GELEI=',N_GELEI,
     &     ', N_HYPRE_ROWS=',N_HYPRE_ROWS,', N_HYPRE_WORK=',N_HYPRE_WORK
         PAUSE
      ENDIF

   13 CONTINUE
      END SUBROUTINE HYPRE_ALLOCATE

C======================================================================
      SUBROUTINE GET_GELEI_LSIZE(NERR)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'hypre.h'
      INCLUDE 'mpif.h'
C
C ASSIGN GELEINDEX(IDIM,JDIM,KDIM) AND LSIZE
C
      INTEGER IERR,NERR,I
      EXTERNAL CALC_LSIZE,CALC_GELEI,CALC_ILOW_IHIGH
!----------------------------------------------------------------------
! Compute local size on this processor across all blocks
!----------------------------------------------------------------------
      LSIZE=0
      CALL CALLWORK(CALC_LSIZE,0)
!----------------------------------------------------------------------
! Distribute all local sizes to all processors
!----------------------------------------------------------------------
        CALL MPI_ALLGATHER(LSIZE,1,MPI_INTEGER,PSIZE,
     &                     1,MPI_INTEGER,
     &                     MPI_COMM_WORLD,IERR)
!----------------------------------------------------------------------
! Reset LSIZE to count sum of local sizes on all previous processors
!----------------------------------------------------------------------
      LSIZE=0
      DO I=0,MYPRC-1
        LSIZE=LSIZE+PSIZE(I)
      ENDDO
!----------------------------------------------------------------------
! Number global element indices in GELEI array
!----------------------------------------------------------------------
      CALL CALLWORK(CALC_GELEI,[1,N_GELEI])
!----------------------------------------------------------------------
! Reset LSIZE to local size, in case it is needed later
!----------------------------------------------------------------------
      LSIZE = PSIZE(MYPRC)
!----------------------------------------------------------------------
! Fill GELEI in ghost layers
!----------------------------------------------------------------------
      CALL UPDATE(N_GELEI,2)
!----------------------------------------------------------------------
! Count lowest and highest row and column indices
!----------------------------------------------------------------------
      ILOW=-1; IHIGH=-1
      CALL CALLWORK(CALC_ILOW_IHIGH,[1,N_GELEI])

      END SUBROUTINE GET_GELEI_LSIZE

C======================================================================
      SUBROUTINE CALC_LSIZE(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'hypre.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),
     &     JL2V(KDIM),KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK
      INTEGER GELEINDEX(IDIM,JDIM,KDIM)
      INTEGER I,J,K

      DO 10 K = 1,KDIM
      DO 10 J = 1,JDIM
      DO 10 I = 1,IDIM
         IF (KEYOUT(I,J,K).NE.1) CYCLE
         LSIZE = LSIZE + 1
   10 CONTINUE

      END SUBROUTINE CALC_LSIZE

C======================================================================
      SUBROUTINE CALC_GELEI(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,GELEI)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'hypre.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),
     &   JL2V(KDIM),KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK,
     &   GELEI(IDIM,JDIM,KDIM)
      INTEGER I,J,K

      GELEI(:,:,:) = 0
      DO 10 K = 1,KDIM
      DO 10 J = 1,JDIM
      DO 10 I = 1,IDIM
         IF (KEYOUT(I,J,K).NE.1) CYCLE
         LSIZE = LSIZE + 1
         GELEI(I,J,K) = LSIZE
   10 CONTINUE

      END SUBROUTINE CALC_GELEI

C======================================================================
      SUBROUTINE CALC_ILOW_IHIGH(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,GELEI)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'hypre.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),
     &   JL2V(KDIM),KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK,
     &   GELEI(IDIM,JDIM,KDIM)
      INTEGER I,J,K

!----------------------------------------------------------------------
! Get highest/lowest row index from active cells
!----------------------------------------------------------------------
      DO 10 K = 1,KDIM
      DO 10 J = 1,JDIM
      DO 10 I = 1,IDIM
        IF (KEYOUT(I,J,K).NE.1) CYCLE
        IF (ILOW.GT.GELEI(I,J,K).OR.ILOW.EQ.-1) ILOW=GELEI(I,J,K)
        IF (IHIGH.LT.GELEI(I,J,K).OR.IHIGH.EQ.-1) IHIGH=GELEI(I,J,K)
   10 CONTINUE

! bag8 debug
      IF (DBGEV) THEN
      WRITE(*,'(2A,2I,3(A,I))')'In CALC_ILOW_IHIGH, ',
     &   'MYPRC,NBLK=',MYPRC,NBLK,
     &   ': LSIZE=',LSIZE,', ILOW=',ILOW,', IHIGH=',IHIGH
      WRITE(10+MYPRC,'(2A,2I,3(A,I))')'In CALC_ILOW_IHIGH, ',
     &   'MYPRC,NBLK=',MYPRC,NBLK,
     &   ': LSIZE=',LSIZE,', ILOW=',ILOW,', IHIGH=',IHIGH
      WRITE(10+MYPRC,*)'GELEI='
      I=IL1
      DO K = KL1-1,KL2+1
      DO J = JL1V(K)-1,JL2V(K)+1
        IF (GELEI(I,J,K).GT.0) THEN
          IF (KEYOUT(I,J,K).EQ.1) THEN
            WRITE(10+MYPRC,'(A,I,A,$)')' ',GELEI(I,J,K),' '
          ELSE
            WRITE(10+MYPRC,'(A,I,A,$)')'(',GELEI(I,J,K),')'
          ENDIF
        ELSE
          WRITE(10+MYPRC,'(A,$)')'     '
        ENDIF
      ENDDO
      WRITE(10+MYPRC,*)
      ENDDO
      ENDIF

      END SUBROUTINE CALC_ILOW_IHIGH

C======================================================================
      SUBROUTINE GETCOLSVALS_7PT(I,J,K,COLS,VALUES,NNZ,KEYOUT,
     &                            GELEINDEX,COF,IDIM,JDIM,KDIM)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
C
C     DIAGONAL ENTRY STARTS FIRST
C
      INTEGER I,J,K,NNZ,IDIM,JDIM,KDIM
      INTEGER COLS(7)
      DOUBLE PRECISION VALUES(7)
      INTEGER KEYOUT(IDIM,JDIM,KDIM),GELEINDEX(IDIM,JDIM,KDIM)
      REAL*4 COF(IDIM,JDIM,KDIM,7)  !!!! DIfferent than in GETCOLSVALS_27PT

      INTEGER OFFSET(3,7),II,JJ,KK,INDEX
      DATA OFFSET/ 0, 0, 0,
     &            -1, 0, 0,  1, 0, 0,    !  X-,  X+
     &             0,-1, 0,  0, 1, 0,    !  Y-,  Y+
     &             0, 0,-1,  0, 0, 1 /   !  Z-,  Z+
      LOGICAL :: DBG = .FALSE.

      COLS = 0
      VALUES = 0.D0
      NNZ = 0

      DO INDEX = 1,7
         II = I + OFFSET(1,INDEX)
         JJ = J + OFFSET(2,INDEX)
         KK = K + OFFSET(3,INDEX)

         IF (KEYOUT(II,JJ,KK).NE.0) THEN
             NNZ = NNZ + 1
             COLS(NNZ) = GELEINDEX(II,JJ,KK)
             VALUES(NNZ) = COF(I,J,K,INDEX)
         ENDIF
      ENDDO

      IF (DBG) THEN
         WRITE(0,'(A,3I3,A,I2)')'In GETCOLSVALS_7PT, I,J,K',I,J,K,
     &                     ', NNZ=',NNZ-1
         WRITE(0,*)'    COF=',COF(I,J,K,1:7)
         WRITE(0,*)'    COLS=',COLS
         WRITE(0,*)'    VALUES=',VALUES
         PAUSE
      ENDIF

      END SUBROUTINE GETCOLSVALS_7PT

C======================================================================
      SUBROUTINE GETCOLSVALS_7PT_R8(I,J,K,COLS,VALUES,NNZ,KEYOUT,
     &                            GELEINDEX,COF,IDIM,JDIM,KDIM)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
C
C     DIAGONAL ENTRY STARTS FIRST
C
      INTEGER I,J,K,NNZ,IDIM,JDIM,KDIM
      INTEGER COLS(7)
      DOUBLE PRECISION VALUES(7)
      INTEGER KEYOUT(IDIM,JDIM,KDIM),GELEINDEX(IDIM,JDIM,KDIM)
      REAL*8 COF(IDIM,JDIM,KDIM,7)  !!!! DIfferent than in GETCOLSVALS_27PT

      INTEGER OFFSET(3,7),II,JJ,KK,INDEX
      DATA OFFSET/ 0, 0, 0,
     &            -1, 0, 0,  1, 0, 0,    !  X-,  X+
     &             0,-1, 0,  0, 1, 0,    !  Y-,  Y+
     &             0, 0,-1,  0, 0, 1 /   !  Z-,  Z+
      LOGICAL :: DBG = .FALSE.

      COLS = 0
      VALUES = 0.D0
      NNZ = 0

      DO INDEX = 1,7
         II = I + OFFSET(1,INDEX)
         JJ = J + OFFSET(2,INDEX)
         KK = K + OFFSET(3,INDEX)

         IF (KEYOUT(II,JJ,KK).NE.0) THEN
             NNZ = NNZ + 1
             COLS(NNZ) = GELEINDEX(II,JJ,KK)
             VALUES(NNZ) = COF(I,J,K,INDEX)
         ENDIF
      ENDDO

      IF (DBG) THEN
         WRITE(0,'(A,3I3,A,I2)')'In GETCOLSVALS_7PT, I,J,K',I,J,K,
     &                     ', NNZ=',NNZ-1
         WRITE(0,*)'    COF=',COF(I,J,K,1:7)
         WRITE(0,*)'    COLS=',COLS
         WRITE(0,*)'    VALUES=',VALUES
         PAUSE
      ENDIF

      END SUBROUTINE GETCOLSVALS_7PT_R8

c======================================================================
      SUBROUTINE GETCOLSVALS_7PT_2VAR(I,J,K,COLS,VALUES,NNZ,KEYOUT,
     &                            GELEINDEX,COF,IDIM,JDIM,KDIM)
C======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
C
C     DIAGONAL ENTRY STARTS FIRST
C
      INTEGER I,J,K,NNZ,IDIM,JDIM,KDIM
      INTEGER COLS(14)
      DOUBLE PRECISION VALUES(14)
      INTEGER KEYOUT(IDIM,JDIM,KDIM),GELEINDEX(IDIM,JDIM,KDIM)
      REAL*4 COF(IDIM,JDIM,KDIM,7,2,2)

      INTEGER OFFSET(3,7),II,JJ,KK,INDEX,L
      DATA OFFSET/ 0, 0, 0,
     &            -1, 0, 0,  1, 0, 0,    !  X-,  X+
     &             0,-1, 0,  0, 1, 0,    !  Y-,  Y+
     &             0, 0,-1,  0, 0, 1 /   !  Z-,  Z+
      LOGICAL :: DBG = .FALSE.

      STOP 'GETCOLSVALS_7PT_2VAR is NOT FINISHED'

      COLS = 0
      VALUES = 0.D0
      NNZ = 0

      DO INDEX = 1,7
         II = I + OFFSET(1,INDEX)
         JJ = J + OFFSET(2,INDEX)
         KK = K + OFFSET(3,INDEX)

         IF (KEYOUT(II,JJ,KK).NE.0) THEN
             NNZ = NNZ + 1
             COLS(NNZ) = GELEINDEX(II,JJ,KK)
ccc             VALUES(NNZ) = COF(I,J,K,INDEX,?,?)
         ENDIF
      ENDDO

      IF (DBG) THEN
         WRITE(0,'(A,3I3,A,I2)')'IN GETCOLSVALS_7PT_2VAR, I,J,K',I,J,K,
     &                     ', NNZ=',NNZ-1
         DO L = 1,7
            WRITE(0,'(A,I3,4E23.15)')'    L,COF=',COF(I,J,K,L,1:2,1:2)
         ENDDO
         WRITE(0,*)'    COLS=',COLS
         WRITE(0,*)'    VALUES=',VALUES
         PAUSE
      ENDIF

      END SUBROUTINE GETCOLSVALS_7PT_2VAR

C======================================================================
C gp: Ideally this routine shoud NOT be sized with MPFA since
C     even from brick we might need it.
C     But since MAPPING is defined in mpfa/fcns.df and extensively used
C     for all MFMFE models, we will keep it there
C======================================================================
C      SUBROUTINE GETCOLSVALS_27PT(I,J,K,COLS,VALUES,NNZ,KEYOUT,
C     &                            GELEI,COF,IDIM,JDIM,KDIM)
C======================================================================
C      IMPLICIT NONE
C      INCLUDE 'control.h'
CC
CC     DIAGONAL ENTRY STARTS FIRST
CC
C      INTEGER I,J,K,NNZ,IDIM,JDIM,KDIM
C      INTEGER COLS(27)
C      DOUBLE PRECISION VALUES(27)
C      INTEGER KEYOUT(IDIM,JDIM,KDIM),GELEI(IDIM,JDIM,KDIM)
C      REAL*8 COF(IDIM,JDIM,KDIM,-13:13)
C
C      INTEGER MAPPING
C
C      INTEGER OFFSET(3,-13:13),II,JJ,KK,INDEX
C      DATA OFFSET/-1,-1,-1,  0,-1,-1,  1,-1,-1,
C     &            -1, 0,-1,  0, 0,-1,  1, 0,-1,
C     &            -1, 1,-1,  0, 1,-1,  1, 1,-1,
C     &            -1,-1, 0,  0,-1, 0,  1,-1, 0,
C     &            -1, 0, 0,  0, 0, 0,  1, 0, 0,
C     &            -1, 1, 0,  0, 1, 0,  1, 1, 0,
C     &            -1,-1, 1,  0,-1, 1,  1,-1, 1,
C     &            -1, 0, 1,  0, 0, 1,  1, 0, 1,
C     &            -1, 1, 1,  0, 1, 1,  1, 1, 1/
C      LOGICAL :: DBG = .FALSE.
C
C      COLS = 0
C      VALUES = 0.D0
C
C      NNZ = 1
C      COLS(NNZ) = GELEI(I,J,K)
C      VALUES(NNZ) = COF(I,J,K,MAPPING(0))
C
C      DO INDEX = -13,13
C         IF (INDEX.EQ.0) CYCLE
C         II = I + OFFSET(1,INDEX)
C         JJ = J + OFFSET(2,INDEX)
C         KK = K + OFFSET(3,INDEX)
C
C         IF (KEYOUT(II,JJ,KK).NE.0) THEN
C             NNZ = NNZ + 1
C             COLS(NNZ) = GELEI(II,JJ,KK)
C             VALUES(NNZ) = COF(I,J,K,MAPPING(INDEX))
C         ENDIF
C      ENDDO
C
C      IF (DBG) THEN
C         WRITE(0,'(A,3I3,A,I2)')'In GETCOLSVALS_27PT, I,J,K',I,J,K,
C     &                     ', NNZ=',NNZ
C         WRITE(0,*)'    COF=',COF(I,J,K,:)
C         WRITE(0,*)'    COLS=',COLS
C         WRITE(0,*)'    VALUES=',VALUES
C         PAUSE
C      ENDIF
C
C      END SUBROUTINE GETCOLSVALS_27PT

c======================================================================
      SUBROUTINE HYPRE_SOLVE(ITLIN,NERR)
c======================================================================
cgp:  Initial version implements flow solve on hexahedra
c     This is now generalized to also handle bricks
c     NOTE: This is a solver for FLOW system;
c     for displacements there is another hypre solver in elastic.f
c======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
C      INCLUDE 'mpfaary.h'
      INCLUDE 'hypre.h'

      INTEGER ITLIN,NERR
      INTEGER COFARR,DUNKARR,RESIDARR,MODACT_OLD

      MODACT_OLD=MODACT
C      MODACT=16
C      MODACT=18
C      MODACT=17

cgp: Added bricks flow models, with 1 variable per cell
C      MODACT=3
cC      MODACT=5
      MODACT=13
C      MODACT=19

      IF (MODACT.EQ.0) THEN
        WRITE(0,*)'hypre.df: MODACT NOT DEFINED'
        STOP
      ENDIF

      COFARR=N_COFV(MODACT)
      DUNKARR=N_DUNKV(MODACT)
      RESIDARR=N_RESIDV(MODACT)

      CALL HYPRE_IPARS_EVFEM(COFARR,DUNKARR,RESIDARR,ITLIN)

      IF (ITLIN>=LSOL_ITMAX) THEN
      NERR = NERR + 1
      IF (MYPRC.EQ.0) THEN
        WRITE(0,'(A,I,A)')'hypre.df: exceeded LSOL_ITMAX=',LSOL_ITMAX,
     &    ' linear iterations'
      ENDIF
      ENDIF

      MODACT = MODACT_OLD

ctm   TAMEEM
C         TOTAL_NUM_LIN_ITER_F = TOTAL_NUM_LIN_ITER_F
C     &            + ITLIN
ctm   TAMEEM

      END SUBROUTINE HYPRE_SOLVE

Cc======================================================================
C      SUBROUTINE HYPRE_TRCHEM(ITLIN,NERR)
Cc======================================================================
Ccgp:  Initial version implements flow solve on hexahedra
Cc     This is now generalized to also handle bricks
Cc     NOTE: This is a solver for FLOW system;
Cc     for displacements there is another hypre solver in elastic.f
C      IMPLICIT NONE
C      INCLUDE 'control.h'
C      INCLUDE 'blkary.h'
CC      INCLUDE 'mpfaary.h'
C      INCLUDE 'hypre.h'
C      INCLUDE 'trarydat.h'
C
C      INTEGER ITLIN,NERR
C      INTEGER COFARR,DUNKARR,RESIDARR
C
C      COFARR=N_TRCOF
C      DUNKARR=N_TRDUNK
C      RESIDARR=N_TRRESID
C      CALL HYPRE_IPARS_EVFEM(COFARR,DUNKARR,RESIDARR,ITLIN)
C
C      IF (ITLIN>=LSOL_ITMAX) THEN
C        NERR = NERR + 1
C        IF (MYPRC.EQ.0) THEN
C          WRITE(0,'(A,I,A)')'hypre.df: exceeded LSOL_ITMAX=',LSOL_ITMAX,
C     &      ' linear iterations'
C        ENDIF
C      ENDIF
C
C      END SUBROUTINE HYPRE_TRCHEM

!======================================================================
! Ben Ganis - Solve the fully coupled EVFEM system.
!======================================================================
      SUBROUTINE HYPRE_IPARS_EVFEM(COFARR,DUNKARR,RESIDARR,ITLIN)
!======================================================================
      IMPLICIT NONE
      INCLUDE 'control.h'
      INCLUDE 'hypre.h'
      INCLUDE 'mpif.h'
! saumik
      INCLUDE 'blkary.h'
! saumik
      INCLUDE 'layout.h'
      INTEGER COFARR,DUNKARR,RESIDARR,ITLIN
      INTEGER IERR,N
      EXTERNAL HYPRE_EVFEM_FILL,HYPRE_EVFEM_SOLN
      double precision res_norm
      integer*8 precond,precond_gotten

      IF (DBGEV)
     &  WRITE(*,'(3(A,I))')
     &    'Begin HYPRE_IPARS_EVFEM, MYPRC=',MYPRC,', ILOW=',ILOW,
     &    ', IHIGH=',IHIGH

      IERR=0
! saumik
      IF(MBPOROE) THEN
! saumik - non-flow processes return without doing anything
         IF(.NOT.MYBLK(2)) RETURN
! saumik - assign flow communicator
         MPI_COMM = FLOW_COMM(1)
      ELSE
         MPI_COMM = MPI_COMM_WORLD
      ENDIF

!----------------------------------------------------------------------
! Create A matrix
!----------------------------------------------------------------------
      call HYPRE_IJMatrixCreate(mpi_comm,ILOW,IHIGH,ILOW,IHIGH,
     &  A,ierr)
      call HYPRE_IJMatrixSetObjectType(A,HYPRE_PARCSR,ierr)
      call HYPRE_IJMatrixInitialize(A,ierr)
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error creating A'
!----------------------------------------------------------------------
! Create b vector
!----------------------------------------------------------------------
      call HYPRE_IJVectorCreate(mpi_comm,ILOW,IHIGH,b,ierr)
      call HYPRE_IJVectorSetObjectType(b,HYPRE_PARCSR,ierr)
      call HYPRE_IJVectorInitialize(b,ierr)
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error creating b'
!----------------------------------------------------------------------
! Create x vector
!----------------------------------------------------------------------
      call HYPRE_IJVectorCreate(mpi_comm,ILOW,IHIGH,x,ierr)
      call HYPRE_IJVectorSetObjectType(x,HYPRE_PARCSR,ierr)
      call HYPRE_IJVectorInitialize(x,ierr)
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error creating x'

!----------------------------------------------------------------------
! Fill A and b
!----------------------------------------------------------------------
      CALL CALLWORK(HYPRE_EVFEM_FILL,[5,COFARR,
     &  RESIDARR,N_GELEI,N_HYPRE_ROWS,N_HYPRE_WORK])
!----------------------------------------------------------------------
! Finish assembly
!----------------------------------------------------------------------
      call HYPRE_IJMatrixAssemble(A,ierr)
      call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
      call HYPRE_IJVectorAssemble(b,ierr)
      call HYPRE_IJVectorGetObject(b,par_b,ierr)
      call HYPRE_IJVectorAssemble(x,ierr)
      call HYPRE_IJVectorGetObject(x,par_x,ierr)
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error assembling'
      IF (DBGEV) call HYPRE_IJMatrixPrint(A,"ij.out.A",ierr)
      IF (DBGEV) call HYPRE_IJVectorPrint(b,"ij.out.b", ierr)
!----------------------------------------------------------------------
! Solve system - AMG
!----------------------------------------------------------------------
      if (hypre_sol_id.eq.0) then
         call HYPRE_BoomerAMGCreate(solver,ierr)
         call HYPRE_BoomerAMGSetPrintLevel(solver,print_level,ierr)
         call HYPRE_BoomerAMGSetCoarsenType(solver,coarsen_type,ierr)
         call HYPRE_BoomerAMGSetRelaxType(solver,relax_type,ierr)
         call HYPRE_BoomerAMGSetMaxIter(solver,lsol_itmax,ierr)
         call HYPRE_BoomerAMGSetNumSweeps(solver,num_sweeps,ierr)
         call HYPRE_BoomerAMGSetMaxLevels(solver,max_levels,ierr)
         call HYPRE_BoomerAMGSetTol(solver,lsol_tol,ierr)
         call HYPRE_BoomerAMGSetStrongThrshld(solver,strong_thres,ierr)
         call HYPRE_BoomerAMGSetCycleType(solver,cycle_type,ierr)
         call HYPRE_BoomerAMGSetup(solver,parcsr_A,par_b,par_x,ierr)
         call HYPRE_BoomerAMGSolve(solver,parcsr_A,par_b,par_x,ierr)
         call HYPRE_BoomerAMGGetNumIterations(solver,ITLIN,ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver,res_norm,ierr)
         call HYPRE_BoomerAMGDestroy(solver,ierr)
!----------------------------------------------------------------------
! Solve system - CG with AMG preconditioner
!----------------------------------------------------------------------
      elseif (hypre_sol_id.eq.1) then
        call HYPRE_ParCSRPCGCreate(mpi_comm,solver,ierr)
        call HYPRE_ParCSRPCGSetMaxIter(solver,lsol_itmax,ierr)
        call HYPRE_ParCSRPCGSetTol(solver,lsol_tol,ierr)
        call HYPRE_ParCSRPCGSetTwoNorm(solver,1,ierr)
        call HYPRE_ParCSRPCGSetRelChange(solver,0,ierr)
        call HYPRE_ParCSRPCGSetPrintLevel(solver,print_level,ierr)
        call HYPRE_BoomerAMGCreate(precond,ierr)
        call HYPRE_BoomerAMGSetCoarsenType(precond,coarsen_type,ierr)
        call HYPRE_BoomerAMGSetRelaxType(precond,relax_type,ierr)
        call HYPRE_BoomerAMGSetMaxIter(precond,1,ierr)
        call HYPRE_BoomerAMGSetCycleType(precond,cycle_type,ierr)
        call HYPRE_BoomerAMGSetMaxLevels(precond,max_levels,ierr)
        call HYPRE_ParCSRPCGSetPrecond(solver,precond_id,precond,ierr)
        call HYPRE_ParCSRPCGGetPrecond(solver,precond_gotten,ierr)
        call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
        call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
        call HYPRE_ParCSRPCGGetNumIterations(solver,ITLIN,ierr)
        call HYPRE_ParCSRPCGGetFinalRelative(solver,res_norm,ierr)
        call HYPRE_BoomerAMGDestroy(precond,ierr)
        call HYPRE_ParCSRPCGDestroy(solver,ierr)
!----------------------------------------------------------------------
! Solve system - GMRES with AMG preconditioner
!----------------------------------------------------------------------
      elseif (hypre_sol_id.eq.2) then
        call HYPRE_ParCSRGMRESCreate(mpi_comm,solver,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 1, ierr=',ierr
        call HYPRE_ParCSRGMRESSetKDim(solver,k_dim,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 2, ierr=',ierr
        call HYPRE_ParCSRGMRESSetMaxIter(solver,lsol_itmax,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 3, ierr=',ierr
        call HYPRE_ParCSRGMRESSetTol(solver,lsol_tol,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 4, ierr=',ierr
        call HYPRE_ParCSRGMRESSetLogging(solver,1,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 5, ierr=',ierr
        call HYPRE_ParCSRGMRESSetPrintLevel(solver,print_level,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 6, ierr=',ierr
        if (precond_id.eq.0) goto 666
        call HYPRE_BoomerAMGCreate(precond,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 7, ierr=',ierr
        call HYPRE_BoomerAMGSetCoarsenType(precond,coarsen_type,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 8, ierr=',ierr
        call HYPRE_BoomerAMGSetRelaxType(precond,relax_type,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 9, ierr=',ierr
        call HYPRE_BoomerAMGSetMeasureType(precond,measure_type,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 10, ierr=',ierr
        call HYPRE_BoomerAMGSetStrongThrshld(precond,strong_thres,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 11, ierr=',ierr
        call HYPRE_BoomerAMGSetTruncFactor(precond,trunc_factor,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 12, ierr=',ierr
        call HYPRE_BoomerAMGSetNumSweeps(precond,1,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 13, ierr=',ierr
        call HYPRE_BoomerAMGSetTol(precond,0.d0,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 14, ierr=',ierr
        call HYPRE_BoomerAMGSetMaxIter(precond,1,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 15, ierr=',ierr
        call HYPRE_BoomerAMGSetCycleType(precond,cycle_type,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 16, ierr=',ierr
        call HYPRE_BoomerAMGSetMaxLevels(precond,max_levels,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 17, ierr=',ierr
        call HYPRE_ParCSRGMRESSetPrecond(solver,precond_id,precond,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 18, ierr=',ierr
        call HYPRE_ParCSRGMRESGetPrecond(solver,precond_gotten,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 19, ierr=',ierr
        if (precond_gotten.ne.precond) then
          print *, 'HYPRE_ParCSRPCGGetPrecond got bad precond'
          stop
        endif
  666   continue
!        write(*,'(a,5i)')'myprc,solver,parcsr_A,par_b,par_x=',
!     &    myprc,solver,parcsr_A,par_b,par_x
        call MPI_BARRIER(mpi_comm,ierr)
        call HYPRE_ParCSRGMRESSetup(solver,parcsr_A,par_b,par_x,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 20, ierr=',ierr
        call HYPRE_ParCSRGMRESSolve(solver,parcsr_A,par_b,par_x,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 21, ierr=',ierr
        call HYPRE_ParCSRGMRESGetNumIteratio(solver,ITLIN,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 22, ierr=',ierr
        call HYPRE_ParCSRGMRESGetFinalRelati(solver,res_norm,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 23, ierr=',ierr
        if (precond_id.ne.0) call HYPRE_BoomerAMGDestroy(precond,ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 24, ierr=',ierr
        call HYPRE_ParCSRGMRESDestroy(solver, ierr)
        if (ierr.ne.0.and.dbgev) write(*,*)'solver 25, ierr=',ierr
      else
        STOP 'Invalid solver id'
      endif
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error solving system'
!----------------------------------------------------------------------
! Get solution x
!----------------------------------------------------------------------
      IF (DBGEV) CALL HYPRE_IJVectorPrint(x,"ij.out.x", ierr)
      CALL CALLWORK(HYPRE_EVFEM_SOLN,[4,DUNKARR,N_GELEI,
     &  N_HYPRE_ROWS,N_HYPRE_WORK])
!----------------------------------------------------------------------
! Deallocate memory
!----------------------------------------------------------------------
      CALL HYPRE_IJMatrixDestroy(A,ierr)
      CALL HYPRE_IJVectorDestroy(b,ierr)
      CALL HYPRE_IJVectorDestroy(x,ierr)
      if (ierr.ne.0.and.dbgev)
     &  write(*,*)'HYPRE_IPARS_EVFEM error deallocating'

      IF (MYPRC.EQ.0) THEN
        IF (IERR.EQ.0) THEN
          WRITE(*,'(1X,1P,A,I6,A,E11.4)')
     &      'HYPRE_EVFEM: ITER=',ITLIN,' LIN RESID=',RES_NORM
        ELSE
          WRITE(*,'(1X,1P,A,I6,A,E11.4,A,I4,A)')
     &      'HYPRE_EVFEM: ITER=',ITLIN,' LIN RESID=',RES_NORM,
     &      ' (IERR=',IERR,')'
        ENDIF
      ENDIF
      CALL HYPRE_CLEARALLERRORS(IERR)

!      IF (DBGEV) THEN
!        WRITE(*,*)'PAUSE at end of HYPRE_IPARS_EVFEM'
!        PAUSE
!      ENDIF

      END SUBROUTINE HYPRE_IPARS_EVFEM

!======================================================================
! Fills the A, b, x HYPRE data structures.
!======================================================================
      SUBROUTINE HYPRE_EVFEM_FILL(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &  JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,COF,RESID,GELEI,ROWS,WORKSOL)
!======================================================================
      IMPLICIT NONE
      include 'control.h'
      include 'layout.h'
C      include 'mpfaary.h'
      include 'hypre.h'
      include 'hypre_dual.h'
      include 'sblkc.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),JL2V(KDIM),
     &     KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK,ITLIN
      INTEGER GELEI(IDIM,JDIM,KDIM)
      REAL*8 RESID(IDIM,JDIM,KDIM),DUNK(IDIM,JDIM,KDIM)
      REAL*4 COF(IDIM,JDIM,KDIM,*)
      INTEGER ROWS(IDIM*JDIM*KDIM)
      DOUBLE PRECISION WORKSOL(IDIM*JDIM*KDIM)

      INTEGER I,J,K,L,IERR,IROW,NNZ,J1,J2,K1,K2,IA,JA,KA
      INTEGER COLS(7)
      DOUBLE PRECISION VALUES(7)

!----------------------------------------------------------------------
! set matrix values from COF on diagonal blocks
!----------------------------------------------------------------------
      DO 10 K=KL1,KL2
      DO 10 J=JL1V(K),JL2V(K)
      DO 10 I=IL1,IL2
         IF(KEYOUT(I,J,K).NE.1) CYCLE

C         IF (MODELON(16)) CALL GETCOLSVALS_27PT(I,J,K,
C     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)
C         IF (MODELON(18)) CALL GETCOLSVALS_27PT(I,J,K,
C     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)
C         IF (MODELON(17)) CALL GETCOLSVALS_27PT(I,J,K,
C     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)

C         IF (MODELON(3)) CALL GETCOLSVALS_7PT(I,J,K,
C     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)
cC         IF (MODELON(5)) CALL GETCOLSVALS_7PT_2VAR(I,J,K,
cC     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)
         IF (MODELON(13)) CALL GETCOLSVALS_7PT(I,J,K,
     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)
C         IF (MODELON(19)) CALL GETCOLSVALS_7PT(I,J,K,
C     &           COLS,VALUES,NNZ,KEYOUT,GELEI,COF,IDIM,JDIM,KDIM)

         IROW = GELEI(I,J,K)
         call HYPRE_IJMatrixSetValues(A,1,NNZ,IROW,COLS,VALUES,IERR)

   10 CONTINUE
      if (ierr.ne.0.and.dbgev) write(*,*)'HYPRE_EVFEM_FILL error 1'

!----------------------------------------------------------------------
! set matrix values from COFINF on off diagonal blocks
!----------------------------------------------------------------------
      IF (HYPRE_EVFEM) THEN
      K1=IIEBS(NBLK)
      K2=K1+NIEBS(NBLK)-1
      DO K=K1,K2         ! Loop over A-block elements K
      IA=IJKS(1,K)       ! Subdomain element on interface
      JA=IJKS(2,K)
      KA=IJKS(3,K)
      J1=ICGES(K)
      J2=J1+NCGES(K)-1
      DO J=J1,J2         ! Loop over A-B interactions J
      IROW=GELEI(IA,JA,KA)
      COLS=0; COLS(1)=GELEI_EV(J,1,1)
      VALUES=0; VALUES(1)=COFINF(J,1,1)
      CALL HYPRE_IJMatrixSetValues(A,1,1,IROW,COLS,VALUES,IERR)
      ENDDO  ! End do-loop over A-B interactions J
      ENDDO  ! End do-loop over A-block faces K
      if (ierr.ne.0.and.dbgev) write(*,*)'HYPRE_EVFEM_FILL error 2'
      ENDIF

!----------------------------------------------------------------------
! set rhs values from RESID
!----------------------------------------------------------------------
      WORKSOL = 0.D0
      ROWS = 0
      IROW = 0
      DO 20 K=KL1,KL2
      DO 20 J=JL1V(K),JL2V(K)
      DO 20 I=IL1,IL2
         IF(KEYOUT(I,J,K).NE.1) CYCLE
         IROW = IROW + 1
         WORKSOL(IROW) = RESID(I,J,K)
         ROWS(IROW) = GELEI(I,J,K)
   20 CONTINUE
      CALL HYPRE_IJVectorSetValues(B, IROW, ROWS, WORKSOL, IERR)
      if (ierr.ne.0.and.dbgev) write(*,*)'HYPRE_EVFEM_FILL error 3'

!----------------------------------------------------------------------
! set solution to be zero
!----------------------------------------------------------------------
      WORKSOL = 0.d0
      CALL HYPRE_IJVectorSetValues(X, IROW, ROWS, WORKSOL, IERR)
      if (ierr.ne.0.and.dbgev) write(*,*)'HYPRE_EVFEM_FILL error 4'

      END SUBROUTINE HYPRE_EVFEM_FILL

!======================================================================
! Stores a piece of the global solution into DUNK.
!======================================================================
      SUBROUTINE HYPRE_EVFEM_SOLN(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &     JL1V,JL2V,KL1,KL2,KEYOUT,NBLK,DUNK,GELEI,ROWS,WORKSOL)
!======================================================================
      IMPLICIT NONE
      include 'control.h'
      include 'layout.h'
      include 'hypre.h'

      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V(KDIM),JL2V(KDIM),
     &     KL1,KL2,KEYOUT(IDIM,JDIM,KDIM),NBLK
      REAL*8 DUNK(IDIM,JDIM,KDIM)
      INTEGER GELEI(IDIM,JDIM,KDIM),ROWS(IDIM*JDIM*KDIM)
      DOUBLE PRECISION WORKSOL(IDIM*JDIM*KDIM)

      INTEGER I,J,K,IROW,IERR

      IROW=0
      DO 10 K=KL1,KL2
      DO 10 J=JL1V(K),JL2V(K)
      DO 10 I=IL1,IL2
        IF (KEYOUT(I,J,K).NE.1) CYCLE
        IROW=IROW+1
        ROWS(IROW) = GELEI(I,J,K)
   10 CONTINUE

      IERR=0
      CALL HYPRE_IJVectorGetValues(x,IROW,ROWS,WORKSOL,IERR)
      IF (IERR.NE.0.and.dbgev) write(*,*)'Problem in HYPRE_EVFEM_SOLN'

      IROW=0
      DO 20 K=KL1,KL2
      DO 20 J=JL1V(K),JL2V(K)
      DO 20 I=IL1,IL2
        IF (KEYOUT(I,J,K).NE.1) CYCLE
        IROW=IROW+1
        DUNK(I,J,K) = WORKSOL(IROW)
   20 CONTINUE

      END SUBROUTINE HYPRE_EVFEM_SOLN
