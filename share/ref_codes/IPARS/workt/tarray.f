C  TARRAY.F - CREATE IMPLICIT SINGLE PHASE FLOW MODEL GRID-ELEMENT ARRAYS

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE TARRAY (KERR)

C  CODE HISTORY:

C  JOHN WHEELER     4/29/97  ALPHA CODE
C  BAHAREH MOMKEN   2/15/99  Hydrology-IMPES garray.df is used as
C                            template
C  JOHN WHEELER    04/03/99  IMPLICIT SINGLE PHASE MODEL
C*********************************************************************
      SUBROUTINE TARRAY (KERR)
C*********************************************************************

C  Creates single phase flow  model grid-element arrays

C  KERR = ERROR NUMBER (OUTPUT, INTEGER)

C  NOTE: See tarydat.h for array descriptions

C*********************************************************************
      IMPLICIT NONE
C        INCLUDE 'msjunk.h'
      INCLUDE 'control.h'
      INCLUDE 'blkary.h'
      INCLUDE 'layout.h'
      INCLUDE 'tarydat.h'
      INCLUDE 'terrcalc.h'

      INTEGER KERR,NARR(3)

cgp dbg
      LOGICAL DBG
      DATA DBG /.FALSE./

C  ALLOCATE GRID-ELEMENT ARRAY SPACE

      KERR=0

      CALL ALCGEA ('FLDEN ',2,0,N_FLDEN,KERR)
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('FLDENN ',2,0,N_FLDENN,KERR)
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('PRES ',2,0,N_PRES,KERR)
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('PRESN ',2,0,N_PRESN,KERR)
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('VEL ',2,3,N_VEL,KERR)
      IF (KERR.GT.0) GO TO 99

! bag8 - use same memory for VEL and VX,VY,VZ
      CALL ALCGA_UNPACK(N_VEL,NARR,KERR)
      N_VX=NARR(1)
      N_VY=NARR(2)
      N_VZ=NARR(3)
      CALL CHG_ANAME(N_VX,'VX ')
      CALL CHG_ANAME(N_VY,'VY ')
      CALL CHG_ANAME(N_VZ,'VZ ')
      IF(KERR.GT.0) GO TO 99

!      CALL ALCGEA ('VX ',2,0,N_VX,KERR)
!      IF(KERR.GT.0) GO TO 99

!      CALL ALCGEA ('VY ',2,0,N_VY,KERR)
!      IF(KERR.GT.0) GO TO 99

!      CALL ALCGEA ('VZ ',2,0,N_VZ,KERR)
!      IF(KERR.GT.0) GO TO 99

C  NOTE THAT THE FRAMEWORK ALSO NEEDS POINTERS TO THE NEXT 3 ARRAYS

      CALL ALCGEA ('COFS ',1,7,N_COF,KERR)
      N_COFV(MODACT)=N_COF
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('RESIDS ',2,1,N_RESID,KERR)
      N_RESIDV(MODACT)=N_RESID
      IF (KERR.GT.0) GO TO 99

      CALL ALCGEA ('DELUNK ',2,1,N_DUNK,KERR)
      N_DUNKV(MODACT)=N_DUNK
      IF (KERR.GT.0) GO TO 99

cgp dbg
      IF (DBG) THEN
         WRITE(0,'(A, 3(A,I4))') 'TARRAY: ','N_COF=',N_COF,
     &               ', N_DUNK=',N_DUNK,', N_RESID=',N_RESID
         PAUSE
      ENDIF
cgp dbg

C bag8 - error against true solution

      IF (ITEST.GT.0) THEN

        CALL ALCGEA ('PRES_ERR ',2,0,N_PRES_ERR,KERR)
        IF (KERR.GT.0) GO TO 99

      ENDIF

! bag8 - adapt grids
      IF ((ADAPT_MODE.EQ.2).OR.(ADAPT_MODE.EQ.3)) THEN
        IF (MYPRC.EQ.0) WRITE(*,'(A,I1,A)')' In TARRAY, ADAPT_MODE=',
     &    ADAPT_MODE,' : allocating extra arrays!'
        CALL ALCGEA ('PRESFINE ',2,0,N_PRESFINE,KERR)
        IF (KERR.GT.0) GO TO 99
      ENDIF

   99 RETURN
      END
