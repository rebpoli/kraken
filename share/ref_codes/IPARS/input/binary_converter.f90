! Utility program to convert ascii to binary and binary to ascii
! created by Ben Ganis 5/24/17.

      PROGRAM binary_converter
      IMPLICIT NONE
      INTEGER :: I,J,K,L,NUMBLKS,IERR,N,ND
      INTEGER, ALLOCATABLE :: NX(:),NY(:),NZ(:),N4(:)
      DOUBLE PRECISION, ALLOCATABLE :: BUF(:,:,:,:)
      CHARACTER*80 :: INFILE, OUTFILE, INPUTNAME
      INTEGER MODE

      WRITE(*,*)'-------------------------------------'
      WRITE(*,*)'Binary input conversion utility      '
      WRITE(*,*)'-------------------------------------'
      WRITE(*,*)'Type 1 to convert ascii to binary'
      WRITE(*,*)'Type 2 to convert binary to ascii'
      WRITE(*,'(a,$)')'Enter mode: '
      READ(*,*) MODE
      IF ((MODE.NE.1).AND.(MODE.NE.2)) STOP 'Error'
      WRITE(*,'(a,$)')'Enter input filename without suffix: '
      READ(*,*) INPUTNAME
      WRITE(*,'(a,$)')'Enter number of fault blocks: '
      READ(*,*) NUMBLKS
      ALLOCATE(NX(NUMBLKS),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate NX'
      ALLOCATE(NY(NUMBLKS),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate NY'
      ALLOCATE(NZ(NUMBLKS),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate NZ'
      ALLOCATE(N4(NUMBLKS),STAT=IERR)
      IF (IERR.NE.0) STOP 'Could not allocate N4'
      DO N=1,NUMBLKS
      WRITE(*,'(a,i2,a,$)')'Enter I,J,K,L dims for block ',N,' : '
      READ(*,*) NX(N),NY(N),NZ(N),N4(N)
      ENDDO

      IF (MODE.EQ.1) THEN

        INFILE=TRIM(INPUTNAME)//'.dat'
        WRITE(*,*)'Reading: ',TRIM(INFILE)
        OPEN(unit=10,file=TRIM(INFILE),status='old')
        OUTFILE=TRIM(INPUTNAME)//'.bin'
        WRITE(*,*)'Writing: ',TRIM(OUTFILE)
        OPEN(unit=11,file=TRIM(OUTFILE),form='binary',status='unknown')

        DO N=1,NUMBLKS
          ALLOCATE(BUF(NX(N),NY(N),NZ(N),N4(N)),STAT=IERR)
          IF (IERR.NE.0) STOP 'Could not allocate BUF'
          READ(10,*) ! comment
          READ(10,*) ((((BUF(I,J,K,L),I=1,NX(N)),J=1,NY(N)),K=1,NZ(N)),&
                         L=1,N4(N))
          DO L=1,N4(N)
          DO K=1,NZ(N)
          DO J=1,NY(N)
          DO I=1,NX(N)
            IF (BUF(I,J,K,L).EQ.0) THEN
              WRITE(*,*)'Error: zero value detected'
              WRITE(*,*)'NBLK,I,J,K,L,VAL=',N,I,J,K,L,BUF(I,J,K,L)
              STOP 1
            ENDIF
          ENDDO
          ENDDO
          ENDDO
          ENDDO

          WRITE(11) BUF
          DEALLOCATE(BUF)
        ENDDO

        CLOSE(10)
        CLOSE(11)

      ELSEIF (MODE.EQ.2) THEN

        INFILE=TRIM(INPUTNAME)//'.bin'
        WRITE(*,*)'Reading: ',TRIM(INFILE)
        OPEN(unit=10,file=TRIM(INFILE),form='binary', &
          access='sequential',status='old')
        OUTFILE=TRIM(INPUTNAME)//'.dat'
        WRITE(*,*)'Writing: ',TRIM(OUTFILE)
        OPEN(unit=11,file=TRIM(OUTFILE), &
             status='unknown')

        DO N=1,NUMBLKS
          ALLOCATE(BUF(NX(N),NY(N),NZ(N),N4(N)),STAT=IERR)
          IF (IERR.NE.0) STOP 'Could not allocate BUF'
          READ(10) BUF
          ND=FLOOR(LOG(REAL(N))/LOG(10.D0))+1  ! Number of digits in N
          WRITE(11,'(2A,I<ND>)')TRIM(INPUTNAME),I,'()='
          WRITE(11,'(1P,6E16.9)') BUF
          DEALLOCATE(BUF)
        ENDDO

        CLOSE(10)
        CLOSE(11)

      ENDIF

      WRITE(*,*)'Done!'
      DEALLOCATE(NX,NY,NZ,N4)

      END PROGRAM binary_converter
