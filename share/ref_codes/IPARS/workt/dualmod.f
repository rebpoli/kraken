!*********************************************************************
!  dualmod.df
!   Ben Ganis 9/22/15  Created dualmod to store nonmatching MFMFE info
!*********************************************************************

      MODULE dualmod
      IMPLICIT NONE
      SAVE

      INTEGER :: NIFNOD,     ! Total number of interface nodes J
     &           NIFACE,     ! Total number of interaction faces K
     &           NORMFACES,  ! Total number of normal faces
     &           NIFNOD1     ! Number of interface nodes before normal nodes
      REAL*8, ALLOCATABLE ::
     &  CIFNOD(:,:),         ! CIFNOD(J,3) = nodal coordinates
     &  DPERMCC(:,:,:),      ! DPERMCC(K,3,3) = 3x3 sym perm tensor, inverse
     &  DPERMINV(:,:,:,:),   ! DPERMINV(J,3,3,8) = inverse perm with density
     &  DPERMINVINI(:,:,:,:),! DPERMINVINI(J,3,3,8) = inverse permeability
     &  DPRES(:,:),          ! DPRES(J,8) = surrounding pressures
     &  DAREARATIO(:)        ! DAREARATIO(K) = ratio of A-B block intf area
                             !                 to A-block element face area
      INTEGER, ALLOCATABLE ::
     &  DFAC2IJK(:,:,:,:),   ! DFAC2IJK(K,L,M,N) = map from K to block IJK node
     &                       !   K = interaction face
     &                       !   L = 1: A-block  2: B-block
     &                       !   M = 1-4:  vertex 1,2,3,4
     &                       !   N = 1-3:  which I,J,K node
     &  DFAC2NOD(:,:),       ! DFAC2NOD(K,4) = 4 intf nodes J for given
     &                       !   interaction face K
     &  DVOLPROP(:,:),       ! DVOLPROP(J,8) = volume properties
     &  DFACEPROP(:,:),      ! DFACEPROP(J,12) = face properties
     &  DFAC2HEX(:,:),       ! DFAC2HEX(K,8) = nodes of B-block hex, face K
     &  DFACENUM(:,:),       ! DFACENUM(J,12) = normal faces at nodes
     &  DNFACE2NOD(:,:),     ! DNFACE2NOD(K,4) = 4 nodes of normal faces
     &  DNOD2FAC(:,:),       ! DNOD2FAC(J,8) = 8 surrounding interaction faces
     &  DLOCALIJK(:,:,:,:)   ! DLOCALIJK(B,J,3,8) = ijk of surrounding 8 vols
                             !                      on block B, node J

      CONTAINS

!*********************************************************************

      SUBROUTINE NEW_IFNOD()
      IMPLICIT NONE
      INTEGER :: IERR
      REAL*8, ALLOCATABLE :: TMP(:,:)

      IF (NIFNOD.GT.0) THEN
        ALLOCATE(TMP(NIFNOD,3),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate TMP in NEW_IFNOD'
        TMP(1:NIFNOD,1:3) = CIFNOD(1:NIFNOD,1:3)
        IF (ALLOCATED(CIFNOD)) DEALLOCATE(CIFNOD)
      ENDIF

      ALLOCATE(CIFNOD(NIFNOD+1,3),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate CIFNOD in NEW_IFNOD'
      IF (NIFNOD.GT.0) THEN
        CIFNOD(1:NIFNOD,1:3)=TMP(1:NIFNOD,1:3)
      ENDIF
      CIFNOD(NIFNOD+1,:)=0

      NIFNOD=NIFNOD+1
      IF (ALLOCATED(TMP)) DEALLOCATE(TMP)

      END SUBROUTINE NEW_IFNOD

!*********************************************************************
! bag8 - add interface node X,Y,Z coordinate to CIFNOD array if not
!        already present, and return the index of the node.
!*********************************************************************
      SUBROUTINE ADD_IFNOD(X,Y,Z,INOD)
!*********************************************************************
      IMPLICIT NONE
      INCLUDE 'sblkc.h'
      REAL*8, INTENT(IN) :: X,Y,Z   ! Input: node coordinate
      INTEGER, INTENT(OUT) :: INOD  ! Output: node index
      INTEGER I
      REAL*8 D
      REAL*8, PARAMETER :: DTOL = 1.D-6

! Check if interface node is already stored (slow...)
      DO I=1,NIFNOD
        D=SQRT((X-CIFNOD(I,1))**2 + (Y-CIFNOD(I,2))**2 +
     &         (Z-CIFNOD(I,3))**2)
        IF (D.LT.DTOL) GOTO 10
      ENDDO

! Node coord is not already present; store in new CIFNOD entry
      CALL NEW_IFNOD()    ! This increments NIFNOD
      CIFNOD(NIFNOD,1:3) = [X,Y,Z]
      INOD = NIFNOD
      RETURN

! Node already present
 10   CONTINUE
      INOD = I
      RETURN

      END SUBROUTINE ADD_IFNOD

!*********************************************************************

      SUBROUTINE NEW_IFACE()
      IMPLICIT NONE
      INTEGER :: IERR
      REAL*8, ALLOCATABLE :: TMP(:,:,:,:),TMP2(:,:),TMP3(:)

      IF (NIFACE.GT.0) THEN
        ALLOCATE(TMP(NIFACE,2,4,3),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate TMP in NEW_IFACE'
        TMP(:,:,:,:) = DFAC2IJK(:,:,:,:)

        ALLOCATE(TMP2(NIFACE,4),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate TMP2 in NEW_IFACE'
        TMP2(:,:) = DFAC2NOD(:,:)

        ALLOCATE(TMP3(NIFACE),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate TMP3 in NEW_IFACE'
        TMP3(:) = DAREARATIO(:)
      ENDIF
      IF (ALLOCATED(DFAC2IJK)) DEALLOCATE(DFAC2IJK)
      IF (ALLOCATED(DFAC2NOD)) DEALLOCATE(DFAC2NOD)
      IF (ALLOCATED(DAREARATIO)) DEALLOCATE(DAREARATIO)

      ALLOCATE(DFAC2IJK(NIFACE+1,2,4,3),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DFAC2IJK in NEW_IFACE'
      ALLOCATE(DFAC2NOD(NIFACE+1,4),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DFAC2NOD in NEW_IFACE'
      ALLOCATE(DAREARATIO(NIFACE+1),STAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*)'NIFACE,IERR=',NIFACE,IERR
        STOP 'Could not allocate DAREARATIO in NEW_IFACE'
      ENDIF
      IF (NIFACE.GT.0) THEN
        DFAC2IJK(1:NIFACE,:,:,:)=TMP(1:NIFACE,:,:,:)
        DFAC2NOD(1:NIFACE,:)=TMP2(1:NIFACE,:)
        DAREARATIO(1:NIFACE)=TMP3(1:NIFACE)
      ENDIF
      DFAC2IJK(NIFACE+1,:,:,:)=0
      DFAC2NOD(NIFACE+1,:)=0
      DAREARATIO(NIFACE+1)=0

      NIFACE=NIFACE+1
      IF (ALLOCATED(TMP)) DEALLOCATE(TMP)
      IF (ALLOCATED(TMP2)) DEALLOCATE(TMP2)
      IF (ALLOCATED(TMP3)) DEALLOCATE(TMP3)

      END SUBROUTINE NEW_IFACE

!*********************************************************************

      SUBROUTINE NEW_NFACE()
      IMPLICIT NONE
      INTEGER :: IERR
      REAL*8, ALLOCATABLE :: TMP(:,:)

      IF (NORMFACES.GT.0) THEN
        ALLOCATE(TMP(NORMFACES,4),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate TMP in NEW_NFACE'
        TMP(:,:) = DNFACE2NOD(:,:)
        IF (ALLOCATED(DNFACE2NOD)) DEALLOCATE(DNFACE2NOD)
      ENDIF

      ALLOCATE(DNFACE2NOD(NORMFACES+1,4),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DNFACE2NOD in NEW_NFACE'
      IF (NORMFACES.GT.0) THEN
        DNFACE2NOD(1:NORMFACES,:)=TMP(1:NORMFACES,:)
      ENDIF
      DNFACE2NOD(NORMFACES+1,:)=0

      NORMFACES=NORMFACES+1
      IF (ALLOCATED(TMP)) DEALLOCATE(TMP)

      END SUBROUTINE NEW_NFACE

!*********************************************************************
! bag8 - add normal face with nodes K1-K4 to DNFACE2NOD array if not
!        already present, and return the index of the normal face.
!*********************************************************************
      SUBROUTINE ADD_NFACE(K1,K2,K3,K4,NFAC)
!*********************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K1,K2,K3,K4  ! Input: node numbers
      INTEGER, INTENT(OUT) :: NFAC        ! Output: normal face index
      INTEGER I,K(4)

! Check if interface node is already stored
      DO I=1,NORMFACES
        K(1:4)=DNFACE2NOD(I,1:4)
        IF ((K1.NE.K(1)).AND.(K1.NE.K(2)).AND.
     &      (K1.NE.K(3)).AND.(K1.NE.K(4))) CYCLE
        IF ((K2.NE.K(1)).AND.(K2.NE.K(2)).AND.
     &      (K2.NE.K(3)).AND.(K2.NE.K(4))) CYCLE
        IF ((K3.NE.K(1)).AND.(K3.NE.K(2)).AND.
     &      (K3.NE.K(3)).AND.(K3.NE.K(4))) CYCLE
        IF ((K4.NE.K(1)).AND.(K4.NE.K(2)).AND.
     &      (K4.NE.K(3)).AND.(K4.NE.K(4))) CYCLE
        GOTO 10
      ENDDO

! Normal face not already present
      CALL NEW_NFACE()    ! This increments NORMFACES
      DNFACE2NOD(NORMFACES,1:4) = [K1,K2,K3,K4]
      NFAC = NORMFACES + NIFACE
      RETURN

! Node already present
 10   CONTINUE
      NFAC = I + NIFACE
      RETURN

      END SUBROUTINE ADD_NFACE

!*********************************************************************

      SUBROUTINE FILL_DNOD2FAC(IDIM,JDIM,KDIM,LDIM,IL1,IL2,
     &   JL1V,JL2V,KL1,KL2,KEYOUT,NBLK)
      IMPLICIT NONE
      INCLUDE 'sblkc.h'
      INTEGER :: IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK
      INTEGER :: JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER :: I,J,K,KD1,KD2,N,IA1,JA1,KA1,IA2,JA2,KA2,K1,K2,J1,J2

! Map from local node (1-4) and face dir (1-6) to vol dir (1-8)
      INTEGER MAP(4,6),HI2LO(6)
      DATA MAP/7,6,2,3, 7,8,4,3, 7,8,5,6,
     &         8,5,1,4, 6,5,1,2, 3,4,1,2/
      DATA HI2LO/4,5,6,1,2,3/

      IF (NIEBS(NBLK).EQ.0) RETURN

      K1=IIEBS(NBLK)
      K2=K1+NIEBS(NBLK)-1
      DO K=K1,K2
      J1=ICGES(K)
      J2=J1+NCGES(K)-1

      IA2=IJKS(1,K)  ! Get IJK inside A block
      JA2=IJKS(2,K)
      KA2=IJKS(3,K)

      DO J=J1,J2
      KD1=KDIRS(J)   ! Direction for getting volume outside A block
      KD2=HI2LO(KD1) ! Direction for getting volume inside A block

      IA1=IA2        ! Get IJK outside A block
      JA1=JA2
      KA1=KA2
      IF (KD1.EQ.1) THEN
        IA1=IA1+1
      ELSEIF (KD1.EQ.2) THEN
        JA1=JA1+1
      ELSEIF (KD1.EQ.3) THEN
        KA1=KA1+1
      ELSEIF (KD1.EQ.4) THEN
        IA1=IA1-1
      ELSEIF (KD1.EQ.5) THEN
        JA1=JA1-1
      ELSEIF (KD1.EQ.6) THEN
        KA1=KA1-1
      ENDIF

      DO I=1,4
      N=DFAC2NOD(J,I)
      DNOD2FAC(N,MAP(I,KD1))=J

      DLOCALIJK(NBLK,N,:,MAP(I,KD1))=[IA1,JA1,KA1]
      DLOCALIJK(NBLK,N,:,MAP(I,KD2))=[IA2,JA2,KA2]

      ENDDO  ! End loop over local node I
      ENDDO  ! End loop over interface J
      ENDDO  ! End loop over element K

      END SUBROUTINE FILL_DNOD2FAC

!*********************************************************************

      SUBROUTINE ALLOC_DUALMOD
      IMPLICIT NONE
      INCLUDE 'layout.h'
      INTEGER :: IERR

      IF (NIFNOD.GT.0) THEN
        ALLOCATE(DPERMINV(NIFNOD,3,3,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DPERMINV'
        ALLOCATE(DPERMINVINI(NIFNOD,3,3,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DPERMINVINI'
        ALLOCATE(DPRES(NIFNOD,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DPRES'
        ALLOCATE(DVOLPROP(NIFNOD,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DVOLPROP'
        DVOLPROP(:,:) = -1
        ALLOCATE(DFACEPROP(NIFNOD,12),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DFACEPROP'
        DFACEPROP(:,:) = -1
        ALLOCATE(DFACENUM(NIFNOD,12),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DFACENUM'
        DFACENUM(:,:) = 0
        ALLOCATE(DFAC2HEX(NIFACE,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DFAC2HEX'
        ALLOCATE(DPERMCC(NIFACE,3,3),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DPERMCC'
        ALLOCATE(DNOD2FAC(NIFNOD,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DNOD2FAC'
        DNOD2FAC(:,:) = 0
        ALLOCATE(DLOCALIJK(NUMBLK,NIFNOD,3,8),STAT=IERR)
        IF (IERR.NE.0) STOP 'Could not allocate DLOCALIJK'
        DLOCALIJK(:,:,:,:) = 1
      ENDIF
      NIFNOD1 = NIFNOD

      END SUBROUTINE ALLOC_DUALMOD

!*********************************************************************

      SUBROUTINE REALLOC_DUALMOD
      IMPLICIT NONE
      INTEGER :: IERR
      INTEGER, ALLOCATABLE :: ITMP(:,:)
      REAL*8, ALLOCATABLE :: DTMP(:,:),DTMP2(:,:,:,:)

      IF (NIFNOD.EQ.0) RETURN

      ALLOCATE(DTMP2(NIFNOD1,3,3,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DTMP2'
      DTMP2=DPERMINV
      DEALLOCATE(DPERMINV)
      ALLOCATE(DPERMINV(NIFNOD,3,3,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DPERMINV'
      DPERMINV(1:NIFNOD1,:,:,:) = DTMP2(:,:,:,:)
      DEALLOCATE(DTMP2)

      ALLOCATE(DTMP2(NIFNOD1,3,3,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DTMP2'
      DTMP2=DPERMINVINI
      DEALLOCATE(DPERMINVINI)
      ALLOCATE(DPERMINVINI(NIFNOD,3,3,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DPERMINVINI'
      DPERMINVINI(1:NIFNOD1,:,:,:) = DTMP2(:,:,:,:)
      DEALLOCATE(DTMP2)

      ALLOCATE(DTMP(NIFNOD1,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DTMP'
      DTMP=DPRES
      DEALLOCATE(DPRES)
      ALLOCATE(DPRES(NIFNOD,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DPRES'
      DPRES(1:NIFNOD1,:) = DTMP(:,:)
      DEALLOCATE(DTMP)

      ALLOCATE(ITMP(NIFNOD1,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate ITMP'
      ITMP=DVOLPROP
      DEALLOCATE(DVOLPROP)
      ALLOCATE(DVOLPROP(NIFNOD,8),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DVOLPROP'
      DVOLPROP(1:NIFNOD1,:) = ITMP(:,:)
      DEALLOCATE(ITMP)

      ALLOCATE(ITMP(NIFNOD1,12),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate ITMP'
      ITMP=DFACEPROP
      DEALLOCATE(DFACEPROP)
      ALLOCATE(DFACEPROP(NIFNOD,12),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DFACEPROP'
      DFACEPROP(1:NIFNOD1,:) = ITMP(:,:)
      DEALLOCATE(ITMP)

      ALLOCATE(ITMP(NIFNOD1,12),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate ITMP'
      ITMP=DFACENUM
      DEALLOCATE(DFACENUM)
      ALLOCATE(DFACENUM(NIFNOD,12),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate DFACENUM'
      DFACENUM(1:NIFNOD1,:) = ITMP(:,:)
      DEALLOCATE(ITMP)

      END SUBROUTINE REALLOC_DUALMOD

!*********************************************************************

      SUBROUTINE DEALLOC_DUALMOD()
      IMPLICIT NONE

      IF (ALLOCATED(CIFNOD))      DEALLOCATE(CIFNOD)
      IF (ALLOCATED(DPERMCC))     DEALLOCATE(DPERMCC)
      IF (ALLOCATED(DPERMINV))    DEALLOCATE(DPERMINV)
      IF (ALLOCATED(DPERMINVINI)) DEALLOCATE(DPERMINVINI)
      IF (ALLOCATED(DPRES))       DEALLOCATE(DPRES)
      IF (ALLOCATED(DFAC2IJK))    DEALLOCATE(DFAC2IJK)
      IF (ALLOCATED(DFAC2NOD))    DEALLOCATE(DFAC2NOD)
      IF (ALLOCATED(DVOLPROP))    DEALLOCATE(DVOLPROP)
      IF (ALLOCATED(DFACEPROP))   DEALLOCATE(DFACEPROP)
      IF (ALLOCATED(DFAC2HEX))    DEALLOCATE(DFAC2HEX)
      IF (ALLOCATED(DFACENUM))    DEALLOCATE(DFACENUM)
      IF (ALLOCATED(DNFACE2NOD))  DEALLOCATE(DNFACE2NOD)
      IF (ALLOCATED(DNOD2FAC))    DEALLOCATE(DNOD2FAC)
      IF (ALLOCATED(DLOCALIJK))   DEALLOCATE(DLOCALIJK)

      END SUBROUTINE DEALLOC_DUALMOD

!*********************************************************************

      END MODULE dualmod

!*********************************************************************

      MODULE debugmod
      IMPLICIT NONE
      SAVE

      INTEGER :: DEBUGMODE = 2,         ! 0=off, 1=ijk, 2=node
     &           DEBUGTRAN = 1,         ! 0=off, 1=on
     &           DEBUGBDRY = 1

      INTEGER, PARAMETER :: NBLK0 = 4,
     &                      IDIM0 = 10,
     &                      JDIM0 = 10,
     &                      KDIM0 = 10
      REAL*8  LCOF_DBG(NBLK0,IDIM0,JDIM0,KDIM0,8,8),
     &        LRHS_DBG(NBLK0,IDIM0,JDIM0,KDIM0,8),
     &        PN_DBG(NBLK0,IDIM0,JDIM0,KDIM0,8),
     &        RHON_DBG(NBLK0,IDIM0,JDIM0,KDIM0,8),
     &        PINV_DBG(NBLK0,IDIM0,JDIM0,KDIM0,3,3,8),
     &        PINVINI_DBG(NBLK0,IDIM0,JDIM0,KDIM0,3,3,8)
      INTEGER NEWT_DBG,
     &        ISTATUS_DBG(NBLK0,IDIM0,JDIM0,KDIM0),
     &        VPROP_DBG(NBLK0,IDIM0,JDIM0,KDIM0,8),
     &        FPROP_DBG(NBLK0,IDIM0,JDIM0,KDIM0,12),
     &        LOCALIJK_DBG(NBLK0,IDIM0,JDIM0,KDIM0,3,8)

      INTEGER, PARAMETER :: DBGNODES = 100
      INTEGER, PARAMETER :: DBGFACES = 100

      REAL*8  LCOF_DBG2(DBGNODES,2,8,8),
     &        LRHS_DBG2(DBGNODES,2,8),
     &        PN_DBG2(DBGNODES,2,8),
     &        RHON_DBG2(DBGNODES,2,8),
     &        PINV_DBG2(DBGNODES,2,3,3,8),
     &        PINVINI_DBG2(DBGNODES,2,3,3,8),
     &        A_DBG2(DBGNODES,2,12,12),
     &        B_DBG2(DBGNODES,2,12,8)
      INTEGER NEWT_DBG2,LASTNOD,LASTHILO,
     &        ISTATUS_DBG2(DBGNODES,2),
     &        VPROP_DBG2(DBGNODES,2,8),
     &        FPROP_DBG2(DBGNODES,2,12),
     &        FDIM_DBG2(DBGNODES,2),
     &        VDIM_DBG2(DBGNODES,2),
     &        LOCALIJK_DBG2(DBGNODES,2,3,8)
      CONTAINS

!*********************************************************************

      SUBROUTINE WRITE_LCOF_DEBUG(E)
      USE dualmod, ONLY : NIFNOD1
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: E
      INTEGER NBLK,I,J,K,L1,L2,L3

      IF (DEBUGMODE.EQ.0) RETURN

! ---------------------- non-matching debug ---------------------------

      IF (DEBUGMODE.EQ.2) THEN
      WRITE(*,'(1X,A,I2)')'Writing LCOF_DEBUG to fort.',20+E

      DO J=1,2
      DO I=1,MIN(DBGNODES,NIFNOD1)

      DO L1=1,8

      DO L2=1,8
        IF (ABS(LCOF_DBG2(I,J,L1,L2)).LT.1.D-10)
     &    LCOF_DBG2(I,J,L1,L2) = 0.D0
      ENDDO

      IF (ABS(LRHS_DBG2(I,J,L1)).LT.1.D-10)
     &  LRHS_DBG2(I,J,L1) = 0.D0
      IF (ABS(PN_DBG2(I,J,L1)).LT.1.D-10)
     &  PN_DBG2(I,J,L1) = 0.D0
      IF (ABS(RHON_DBG2(I,J,L1)).LT.1.D-10)
     &  RHON_DBG2(I,J,L1) = 0.D0

      DO L2=1,3
      DO L3=1,3
        IF (ABS(PINV_DBG2(I,J,L2,L3,L1)).LT.1.D-10)
     &    PINV_DBG2(I,J,L2,L3,L1) = 0.D0
      ENDDO
      ENDDO

      DO L2=1,3
      DO L3=1,3
        IF (ABS(PINVINI_DBG2(I,J,L2,L3,L1)).LT.1.D-10)
     &    PINVINI_DBG2(I,J,L2,L3,L1) = 0.D0
      ENDDO
      ENDDO

      ENDDO

      WRITE(20+E,'(3(A,I3))')'NEWT=',NEWT_DBG2,', NODE=',I,', HILO=',J
      WRITE(20+E,'(A,I1)')'ISTATUS=',ISTATUS_DBG2(I,J)
      WRITE(20+E,'(A,8I3)')'VPROP=',VPROP_DBG2(I,J,:)

!      WRITE(20+E,'(A,12I3)')'FPROP=',FPROP_DBG2(I,J,:)

      WRITE(20+E,'(A,$)')'FPROP='
      DO L1=1,12
        IF ((FPROP_DBG2(I,J,L1).EQ.-1).OR.
     &      (FPROP_DBG2(I,J,L1).EQ.3)) THEN
          WRITE(20+E,'(A,$)')' -1'
        ELSEIF ((FPROP_DBG2(I,J,L1).EQ.6).OR.
     &      (FPROP_DBG2(I,J,L1).EQ.0)) THEN
          WRITE(20+E,'(A,$)')'  6'
        ELSEIF (FPROP_DBG2(I,J,L1).EQ.2) THEN
          WRITE(20+E,'(A,$)')'  2'
        ELSE
          WRITE(20+E,'(I3,$)')FPROP_DBG2(I,J,L1)
        ENDIF
      ENDDO
      WRITE(20+E,*)
      WRITE(20+E,'(A,2I3)')'VDIM,FDIM=',VDIM_DBG2(I,J),FDIM_DBG2(I,J)
      WRITE(20+E,'(A,$)')'LOCALIJK='
      DO L1=1,8
        IF (VPROP_DBG2(I,J,L1).EQ.5) THEN
          WRITE(20+E,'(3I3,$)')LOCALIJK_DBG2(I,J,:,L1)
        ELSE
          WRITE(20+E,'(A,$)')'  -  -  -'
        ENDIF
        IF (L1.LT.8) WRITE(20+E,'(A,$)')','
      ENDDO
      WRITE(20+E,*)
      WRITE(20+E,'(A,8F9.2)')'PN=',PN_DBG2(I,J,:)
      WRITE(20+E,'(A,8F9.2)')'RHON=',RHON_DBG2(I,J,:)
!      WRITE(20+E,'(A)')'PINVINI='
!      WRITE(20+E,'(9F9.2)')PINVINI_DBG2(I,J,:,:,:)
      WRITE(20+E,'(A)')'PINV='
      WRITE(20+E,'(9F9.2)')PINV_DBG2(I,J,:,:,:)
      WRITE(20+E,'(A)')'A='
      DO L1=1,FDIM_DBG2(I,J)
        DO L2=1,FDIM_DBG2(I,J)
          WRITE(20+E,'(F9.2,$)')A_DBG2(I,J,L1,L2)
        ENDDO
        WRITE(20+E,*)
      ENDDO
      WRITE(20+E,'(A)')'B='
      DO L1=1,FDIM_DBG2(I,J)
        DO L2=1,VDIM_DBG2(I,J)
          WRITE(20+E,'(F9.2,$)')B_DBG2(I,J,L1,L2)
        ENDDO
        WRITE(20+E,*)
      ENDDO
      WRITE(20+E,'(A,8F9.2)')'LRHS=',LRHS_DBG2(I,J,:)
      WRITE(20+E,'(A)')'LCOF='
      WRITE(20+E,'(8F9.2)')LCOF_DBG2(I,J,:,:)
      WRITE(20+E,*)

      ENDDO
      ENDDO
!      IF (NEWT_DBG2.EQ.3) STOP 222

! ---------------------- matching debug -------------------------------

      ELSEIF (DEBUGMODE.EQ.1) THEN
      WRITE(*,'(1X,A,I2)')'Writing LCOF_DEBUG to fort.',10+E

      DO NBLK=1,NBLK0
!      DO NBLK=1,1
      DO K=2,KDIM0
      DO J=2,JDIM0
      DO I=1,IDIM0

      DO L1=1,8

      DO L2=1,8
        IF (ABS(LCOF_DBG(NBLK,I,J,K,L1,L2)).LT.1.D-10)
     &    LCOF_DBG(NBLK,I,J,K,L1,L2) = 0.D0
      ENDDO

      IF (ABS(LRHS_DBG(NBLK,I,J,K,L1)).LT.1.D-10)
     &  LRHS_DBG(NBLK,I,J,K,L1) = 0.D0
      IF (ABS(PN_DBG(NBLK,I,J,K,L1)).LT.1.D-10)
     &  PN_DBG(NBLK,I,J,K,L1) = 0.D0
      IF (ABS(RHON_DBG(NBLK,I,J,K,L1)).LT.1.D-10)
     &  RHON_DBG(NBLK,I,J,K,L1) = 0.D0

      DO L2=1,3
      DO L3=1,3
        IF (ABS(PINV_DBG(NBLK,I,J,K,L2,L3,L1)).LT.1.D-10)
     &    PINV_DBG(NBLK,I,J,K,L2,L3,L1) = 0.D0
        IF (ABS(PINVINI_DBG(NBLK,I,J,K,L2,L3,L1)).LT.1.D-10)
     &    PINVINI_DBG(NBLK,I,J,K,L2,L3,L1) = 0.D0
      ENDDO
      ENDDO

      ENDDO

      WRITE(10+E,'(2(A,I3),A,3I3)')'NEWT=',NEWT_DBG,', NBLK=',NBLK,
     &  ', I,J,K=',I,J,K
      WRITE(10+E,'(A,I1)')'ISTATUS=',ISTATUS_DBG(NBLK,I,J,K)
      WRITE(10+E,'(A,8I3)')'VPROP=',VPROP_DBG(NBLK,I,J,K,:)

!      WRITE(10+E,'(A,12I3)')'FPROP=',FPROP_DBG(NBLK,I,J,K,:)

      WRITE(10+E,'(A,$)')'FPROP='
      DO L1=1,12
        IF ((FPROP_DBG(NBLK,I,J,K,L1).EQ.-1).OR.
     &      (FPROP_DBG(NBLK,I,J,K,L1).EQ.3)) THEN
          WRITE(10+E,'(A,$)')' -1'
        ELSEIF ((FPROP_DBG(NBLK,I,J,K,L1).EQ.6).OR.
     &      (FPROP_DBG(NBLK,I,J,K,L1).EQ.0)) THEN
          WRITE(10+E,'(A,$)')'  6'
        ELSEIF (FPROP_DBG(NBLK,I,J,K,L1).EQ.2) THEN
          WRITE(10+E,'(A,$)')'  2'
        ELSE
          WRITE(10+E,'(I3,$)')FPROP_DBG(NBLK,I,J,K,L1)
        ENDIF
      ENDDO
      WRITE(10+E,*)
      WRITE(10+E,'(A,$)')'LOCALIJK='
      DO L1=1,8
        IF (VPROP_DBG(NBLK,I,J,K,L1).EQ.5) THEN
          WRITE(10+E,'(3I3,$)')LOCALIJK_DBG(NBLK,I,J,K,:,L1)
        ELSE
          WRITE(10+E,'(A,$)')'         '
        ENDIF
        IF (L1.LT.8) WRITE(10+E,'(A,$)')','
      ENDDO
      WRITE(10+E,*)
!      WRITE(10+E,'(A,8F9.2)')'PN=',PN_DBG(NBLK,I,J,K,:)
!      WRITE(10+E,'(A,8F9.2)')'RHON=',RHON_DBG(NBLK,I,J,K,:)
!!      WRITE(10+E,'(A)')'PINVINI='
!!      WRITE(10+E,'(9F9.2)')PINVINI_DBG(NBLK,I,J,K,:,:,:)
!      WRITE(10+E,'(A)')'PINV='
!      WRITE(10+E,'(9F9.2)')PINV_DBG(NBLK,I,J,K,:,:,:)
!      WRITE(10+E,'(A,8F9.2)')'LRHS=',LRHS_DBG(NBLK,I,J,K,:)
!      WRITE(10+E,'(A)')'LCOF='
!      WRITE(10+E,'(8F9.2)')LCOF_DBG(NBLK,I,J,K,:,:)
!      WRITE(10+E,*)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!      IF (NEWT_DBG.EQ.3) STOP 222

      ELSE
        STOP 'Unknown DEBUGMODE'
      ENDIF

      END SUBROUTINE WRITE_LCOF_DEBUG

!*********************************************************************

      END MODULE debugmod

