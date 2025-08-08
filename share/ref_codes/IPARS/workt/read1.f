C  READ1.FOR - FREEFORM KEYWORD INPUT PACKAGE

C  ROUTINES IN THIS MODULE:

C  SUBROUTINE READER  (ENDKEY,INPFIL,INCFIL,KODRET)
C  SUBROUTINE PUTERR  (NERR,EMSG,DAT,LDAT,LPOINT)
C  SUBROUTINE GETNUM  (R8,KEY,L,K)
C  SUBROUTINE UNDEF   (NERR)
C  SUBROUTINE GETGRDA (VNAM,KARY,NUMRET,NERR)
C  SUBROUTINE GRDAIN  (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
C                     KEYOUT,NBLK,VAL)
C  SUBROUTINE DEFAULT (DEFUNT)
C  SUBROUTINE GETBLK  (VNAM,BUF,MAXBUF,NDIM,LOC,LEN,NERR)
C  SUBROUTINE SETNBG  ()
C  SUBROUTINE TBLKOUT (BLK,LB)

C  HISTORY:

C  JOHN WHEELER     10/20/95    ORIGINAL BETA CODE
C  JOHN WHEELER      8/27/99    ADD BLOCK PRINT ROUTINE
C  RICK DEAN        11/26/01    INITIALIZE LQ TO FALSE

C  NOTES:

C     1) ERROR NUMBERS 101 TO 130 ARE RESERVED FOR KEY-WORD INPUT

C     2) THE MAXIMUM SIZE OF THE INPUT FILE (EXCLUDING COMMENTS AND
C        UNNECESSARY BLANKS) IS THE PARAMETER MAXCHR

C     3) THE MAXIMUM RECORD LENGTH IS THE PARAMETER RECLEN

C     4) THE /SCRAT1/ COMMON MEMORY CAN BE USED FOR OTHER PURPOSES
C        BEFORE AND AFTER INPUT PROCESSING

C     5) THIS VERSION OF READER CAN EXTRACT A GRID-ELEMENT ARRAY FROM A
C        GLOBAL ARRAY.

C*********************************************************************
      SUBROUTINE READER (ENDKEY,INPFIL,INCFIL,KODRET)
C*********************************************************************
C  ROUTINE READS THE DATA UP TO A TERMINATOR INTO A SINGLE CHARACTER
C  STRING FOR LATER ANALYSIS.

C  ENDKEY = CHARACTER STRING THAT TERMINATES INPUT.  DATA FOLLOWING
C           THE TERMINATOR IN THE SAME RECORD WILL BE DISCARDED.
C           (INPUT, CHARACTER*20)

C  INPFIL = INPUT FILE NUMBER (MUST BE OPEN). (INPUT, INTEGER)

C  INCFIL = INCLUDE FILE NUMBER (MUST BE CLOSED). (INPUT, INTEGER)
C           MUST NOT EQUAL INPFIL

C  KODRET = RETURN CODE. (OUTPUT, INTEGER)
C         = 0 ==> NO ERRORS ENCOUNTERED
C         > 0 ==> NUMBER OF ERRORS ENCOUNTERED

C  NOTES:  A COLON WILL BE PLACED BEFORE EACH VARIABLE NAME
C          UNNECESSARY BLANKS WILL BE DISCARDED
C          COMMENTS WILL BE DISCARDED
C*********************************************************************
      USE scrat1mod

      IMPLICIT NONE
      INTEGER MAXCHR,NRCLEN,KODRET,INPFIL,NFIL,I,J,K,LKEY,LNBC
     &   ,NET,NETB,KK,K1,INCFIL,LTL1,LTLD
      PARAMETER (MAXCHR=10000000,NRCLEN=202)

      LOGICAL LQ
      CHARACTER*1 C(NRCLEN),ENDKEY(20),BLANK,QUOTE,RECEND,TEXEND,AMBER,
     & DOLLAR,COLON,LEFT,RIGHT,COMMA,AAA,ZZZ,CC,TRU(4),FAL(5),BLOCK(5),
     & LBRAC,RBRAC,INCLUD(7),E1(50)
      CHARACTER*8 ENDBK,EB8
      CHARACTER*50 E
      INCLUDE 'control.h'

      EQUIVALENCE (C,EB8),(E1(1),E)

      DATA BLANK/' '/,QUOTE/'"'/,DOLLAR/'$'/,COLON/':'/,LEFT/'('/,
     & RIGHT/')'/,COMMA/','/,AAA/'A'/,ZZZ/'Z'/,TRU/'T','R','U','E'/,
     & FAL/'F','A','L','S','E'/,AMBER/'&'/,LBRAC/'['/,RBRAC/']'/,
     & BLOCK/'B','l','o','c','k'/,ENDBK/'EndBlock'/,
     & INCLUD/'I','n','c','l','u','d','e'/

      KODRET=0
      RECEND=CHAR(30)
      TEXEND=CHAR(31)
      NFIL=INPFIL
      LQ = .FALSE.

C  COUNT CHARACTERS IN THE TERMINATOR

      DO 1 I=20,1,-1
      LKEY=I
      IF (ENDKEY(I).NE.BLANK) GO TO 2
    1 CONTINUE

C  START READ LOOP
C  READ A RECORD

    2 LAST=0
    3 LAST=LAST+1
      A(LAST)=BLANK
      IF (LAST.EQ.MAXCHR) GO TO 13
   12 READ(NFIL,4,END=20) C
    4 FORMAT(202A1)

C  FIND LAST NONBLANK CHARACTER IN THE RECORD

      DO 5 I=NRCLEN,1,-1
      IF (C(I).NE.BLANK) THEN
         LNBC=I
         GO TO 6
      ENDIF
    5 CONTINUE
      GO TO 7

C  LOOK FOR QUOTED STRINGS

    6 LQ=.FALSE.
      NET=0
      NETB=0

      I=0
   18 I=I+1
      IF (I.GT.LNBC) GO TO 7
      CC=C(I)

      IF (CC.EQ.QUOTE) THEN
         LQ=.NOT.LQ
         GO TO 9
      ENDIF
      IF (LQ) GO TO 9

C  LOOK FOR TERMINATOR

      K=I
      DO 10 J=1,LKEY
      IF (K.GT.LNBC.OR.C(K).NE.ENDKEY(J)) GO TO 31
   10 K=K+1
      GO TO 20

C  PROCESS INCLUDE STATEMENTS

   31 K=I
      DO 30 J=1,7
      IF (K.GT.LNBC.OR.C(K).NE.INCLUD(J)) GO TO 11
   30 K=K+1

      KK=K
      DO 32 J=KK,LNBC
      K=J
      IF (C(J).NE.BLANK.AND.C(J).NE.'=') GO TO 34
   32 CONTINUE
      GO TO 33

   34 E=' '
      DO 35 J=1,50
      IF (K.GT.LNBC.OR.C(K).EQ.BLANK) GO TO 36
      E1(J)=C(K)
   35 K=K+1

   36 OPEN (INCFIL,FILE=E,STATUS='OLD',ERR=33)
      NFIL=INCFIL
      GO TO 12

C  LOOK FOR COMMENTS, PARENTHESES, AND BRACKETS

   11 IF (CC.EQ.DOLLAR) GO TO 7
      IF (CC.EQ.LEFT) NET=NET+1
      IF (CC.EQ.RIGHT) NET=NET-1
      IF (CC.EQ.LBRAC) NETB=NETB+1
      IF (CC.EQ.RBRAC) NETB=NETB-1

C  LOOK FOR VARIABLE NAMES (MARK BY A COLON) AND KEY WORDS

      IF (CC.GE.AAA.AND.CC.LE.ZZZ) THEN
C  LOOK FOR BLOCK
         DO 8 J=1,5
         IF (C(J+I-1).NE.BLOCK(J)) GO TO 19
         IF (LAST+J.GE.MAXCHR) GO TO 13
    8    A(LAST+J)=BLOCK(J)
         LAST=LAST+6
         A(LAST)=RECEND
         GO TO 21
C  LOOK FOR TRUE
   19    DO 14 J=1,4
         IF (C(J+I-1).NE.TRU(J)) GO TO 15
         IF (LAST+J.GE.MAXCHR) GO TO 13
   14    A(LAST+J)=TRU(J)
         I=I+3
         LAST=LAST+4
         GO TO 18
C  LOOK FOR FALSE
   15    DO 16 J=1,5
         IF (C(J+I-1).NE.FAL(J)) GO TO 17
         IF (LAST+J.GE.MAXCHR) GO TO 13
   16    A(LAST+J)=FAL(J)
         I=I+4
         LAST=LAST+5
         GO TO 18
C
   17    IF (NET.EQ.0.AND.NETB.EQ.0.AND.(A(LAST).EQ.BLANK.OR.A(LAST)
     &     .EQ.COMMA)) A(LAST)=COLON
      ELSE

C  DISCARD UNNEEDED BLANKS

         IF (CC.EQ.BLANK.AND.(A(LAST).EQ.BLANK.OR.A(LAST).EQ.COMMA))
     &      GO TO 18
         IF (CC.EQ.COMMA.AND.A(LAST).EQ.BLANK) LAST=LAST-1
      ENDIF

    9 LAST=LAST+1
      IF (LAST.GE.MAXCHR) GO TO 13
      A(LAST)=CC

      GO TO 18

C  READ BLOCK TEXT

   21 READ (NFIL,4,END=13) C
C  TEST FOR END OF BLOCK
      IF (EB8.EQ.ENDBK) THEN
         LAST=LAST+1
         A(LAST)=TEXEND
         GO TO 3
      ENDIF
C  FIND LAST CHARACTER, DISCARDING ANY COMMENT
      K=0
      DO 22 J=1,NRCLEN
      IF (C(J).EQ.DOLLAR) GO TO 23
      IF (C(J).NE.BLANK) K=J
   22 CONTINUE
C  TEST FOR EXCESS DATA
   23 IF (LAST+K.GE.MAXCHR) GO TO 25
C  TEST CONTINUATION
      IF (C(1).EQ.AMBER) THEN
         LAST=LAST-1
         K1=2
      ELSE
         K1=1
      ENDIF
C  COPY RECORD
      LAST=LAST+1
      DO 24 J=K1,K
      A(LAST)=C(J)
   24 LAST=LAST+1
      A(LAST)=RECEND
      GO TO 21

C  ERROR MESSAGES

    7 IF (LQ) THEN
         E='UNMATCHED " IN A RECORD'
         CALL PUTERR(101,E,C,LNBC,LNBC)
         A(LAST)=QUOTE
         KODRET=KODRET+1
      ENDIF
      GO TO 3

   33 E='INVALID INCLUDE FILE'
      CALL PUTERR(101,E,C,LNBC,KK)
      KODRET=KODRET+1
      CLOSE (INCFIL)
      GO TO 3

   25 E='BLOCK DATA INPUT NOT TERMINATED'
      CALL PUTERR(119,E,C,K,0)
      KODRET=KODRET+1
      GO TO 20

   13 E='INPUT DATA FILE TOO LONG, LAST RECORD READ WAS'
      CALL PUTERR(102,E,C,LNBC,0)
      KODRET=KODRET+1

   20 IF (NFIL.EQ.INCFIL) THEN
         CLOSE (INCFIL)
         NFIL=INPFIL

         LTL1=MAX(2,LAST-128)
         DO 40 LTLD=LAST,LTL1,-1
         IF (A(LTLD).LT.RECEND) A(LTLD)=BLANK
   40    CONTINUE

         GO TO 12
      ENDIF
      A(LAST)=COLON
      RETURN
      END
C*********************************************************************
      SUBROUTINE PUTERR (NERR,EMSG,DAT,LDAT,LPOINT)
C*********************************************************************
C  ROUTINE DISPLAYS ERROR MESSAGES FOR KEY-WORD INPUT

C  NERR   = ERROR NUMBER (INPUT, INTEGER)

C  EMSG   = ERROR MESSAGE (INPUT, CHARACTER*50)

C  DAT()  = INPUT DATA CONTAINING ERROR (INPUT, CHARACTER*1)

C  LDAT   = LENGTH OF DAT TO BE PRINTED (INPUT, INTEGER)
C           MAY BE TRUNCATED BY : OR CHAR(30) BEFORE LDAT IS REACHED

C  LPOINT = OFFSET IN DAT AT WHICH ERROR WAS DETECTED (INPUT, INTEGER)
C         = 0 ==> DO NOT DISPLAY MARKER
C           NO MARKER IF LPOINT > LDAT

C  NOTES:
C     1)  THIS ROUTINE SETS LEVERR = 2

C     2)  THIS ROUTINE PRINTS ONLY ON PROCESSOR 0 UNLESS LEVELA = .TRUE.

C*********************************************************************
      IMPLICIT NONE
      INTEGER MAXCHR, NRCLEN,NERR,LL,I,L,LDAT,LPOINT
      PARAMETER (MAXCHR=10000000,NRCLEN=202)

      CHARACTER*1 DOL,BLK,COLON,RECEND,LIN(79),DAT(*)
      CHARACTER*50 EMSG

      INCLUDE 'control.h'

      DATA DOL/'$'/,BLK/' '/,COLON/':'/

      LEVERR=2

      IF ((MYPRC.NE.0).AND.(.NOT. LEVELA)) RETURN

      IF ((NUMPRC.GT.1).AND.LEVELC) THEN
         WRITE (NFOUT,6) NERR,MYPRC,EMSG
    6    FORMAT(/' ERROR #',I4,' ON PROCESSOR',I4/1X,A50)
      ELSE
         WRITE (NFOUT,1) NERR,EMSG
    1    FORMAT (/' ERROR #',I4,'; ',A50)
      ENDIF

      LL=0
      RECEND=CHAR(30)
      DO 4 I=1,LDAT
      IF (DAT(I).EQ.COLON.OR.DAT(I).EQ.RECEND) GO TO 5
    4 LL=LDAT
    5 IF (LL.GT.0) THEN
         WRITE (NFOUT,2) (DAT(I),I=1,LL)
    2    FORMAT(1X,79A1)
         IF (LPOINT.NE.0.AND.LPOINT.LE.LL) THEN
            L=MOD(LPOINT,NRCLEN)
            DO 3 I=1,L
    3       LIN(I)=BLK
            LIN(L)=DOL
            WRITE(NFOUT,2) (LIN(I),I=1,L)
         ENDIF
      ENDIF
      RETURN
      END
C*********************************************************************
      SUBROUTINE GETNUM (R8,KEY,L,K)
C*********************************************************************
C   ASCII TO REAL*8 TRANSLATION ROUTINE

C   R8  = NUMBER OUTPUT (OUTPUT, REAL*8)

C   KEY = RESULT KEY (OUTPUT, INTEGER)
C       = 0 ==> NUMBER FOUND
C       = 1 ==> NOT A NUMBER
C       = 2 ==> UNITS ERROR (MESSAGE GENERATED)

C   L   = 1ST OFFSET IN A() TO BE TESTED (INPUT, INTEGER)

C   K   = 1ST OFFSET IN A() FOLLOWING THE NUMBER (OUTPUT, INTEGER)

C   NOTES:

C      1) IF KEY > 0 THEN R8 AND K ARE UNCHANGED

C*********************************************************************
      USE scrat1mod
      PARAMETER (MAXCHR=10000000)

      REAL*8 R8,D,V
      LOGICAL NONUM,PNT,SGN
      CHARACTER*1 NUM(10),BLANK,COMMA,PLUS,MINUS,POINT,EEE,LBRAC,
     & RBRAC,ASTR,COLON
      CHARACTER*50 E

      INCLUDE 'control.h'
      INCLUDE 'readdat.h'
!      INCLUDE 'scrat1.h'

      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/,
     & BLANK/' '/,COMMA/','/,PLUS/'+'/,MINUS/'-'/,POINT/'.'/,EEE/'E'/
     & LBRAC/'['/,RBRAC/']'/,ASTR/'*'/,COLON/':'/

      KEY=1
      J=L

C  SKIP ANY LEADING BLANKS OR COMMAS

      IF (A(J).EQ.BLANK) J=J+1
      IF (A(J).EQ.COMMA) J=J+1

C  GET DECIMAL SIGN

    2 SGN=.FALSE.
      IF (A(J).EQ.PLUS) THEN
         J=J+1
      ELSE
         IF (A(J).EQ.MINUS) THEN
            SGN=.TRUE.
            J=J+1
         ENDIF
      ENDIF

C  GET DECIMAL FRACTION

      V=0.D0
      NONUM=.TRUE.
      PNT=.FALSE.
    3 IF (A(J).EQ.POINT) THEN
         IF (PNT) RETURN
         PNT=.TRUE.
         D=.1D0
         J=J+1
      ENDIF
      DO 4 I=1,10
      IF (A(J).EQ.NUM(I)) THEN
         IF (PNT) THEN
            V=V+D*(I-1)
            D=.1D0*D
         ELSE
            V=10.D0*V+I-1
         ENDIF
         J=J+1
         NONUM=.FALSE.
         GO TO 3
      ENDIF
    4 CONTINUE
      IF (NONUM) RETURN
      IF (SGN) V=-V

C     GET EXPONENT

      IF (A(J).NE.EEE) GO TO 7
      J=J+1
      SGN=.FALSE.
      IF (A(J).EQ.PLUS) THEN
         J=J+1
      ELSE
         IF (A(J).EQ.MINUS) THEN
            SGN=.TRUE.
            J=J+1
         ENDIF
      ENDIF
      N=0
      NONUM=.TRUE.
    5 DO 6 I=1,10
      IF (A(J).EQ.NUM(I)) THEN
         N=10*N+I-1
         J=J+1
         NONUM=.FALSE.
         GO TO 5
      ENDIF
    6 CONTINUE
      IF (NONUM) RETURN
      IF (SGN) N=-N
      V=V*(10.D0**N)

C  GET UNITS

    7 IF (ISUNTD.AND.NOTINDX) THEN
         ISUNT=.TRUE.
         ISUNTD=.FALSE.
         CALL CONVRT(UNTDEF,UNTSTD,FACMU,FACAU,KODRET,E)
         IF (KODRET.NE.0) THEN
            NERR=NERR+1
            ISUNT=.FALSE.
            LEVERR=2
            IF (LEVELC) THEN
               WRITE (NFOUT,14) KODRET,E
               WRITE (NFOUT,15) (UNTDEF(I),I=1,60)
            ENDIF
   14       FORMAT(/' ERROR #',I4,'; ',A50)
   15       FORMAT(' IN DEFAULT UNITS ',60A1)
         ENDIF
      ENDIF

      IF (A(J).EQ.LBRAC) THEN
         LL=LAST-J+1
         IF (LL.GT.60) LL=60
         DO 8 I=1,LL
         JJ=J+I
         IF (A(JJ).EQ.RBRAC) GO TO 9
    8    CONTINUE
         E='MISSING ] IN UNITS SPECFICATION'
         KODRET =121
         GO TO 10
    9    CALL CONVRT(A(J),UNTSTD,FACMU,FACAU,KODRET,E)
         IF (KODRET.NE.0) GO TO 10
         J=JJ+1
         ISUNT=.TRUE.
      ENDIF

C  APPLY UNITS CONVERSION

      IF (ISUNT.AND.NOTINDX.AND.(J.EQ.LAST.OR.A(J+1).NE.ASTR))
     &   V=V*FACMU+FACAU

C  RETURN NUMBER

      KEY=0
      K=J
      R8=V
      RETURN

C  UNITS ERROR

   10 L=1
      LL=LAST-J+1
      IF (LL.GT.65) LL=65
      DO 11 I=1,LL
      M=I
      IF (A(J+I).EQ.COLON) GO TO 13
   11 CONTINUE
   13 CALL PUTERR(KODRET,E,A(J),M,1)
      KEY=2

      END
C*********************************************************************
      SUBROUTINE UNDEF (NERR)
C*********************************************************************
C  ROUTINE CHECKS INPUT DATA FOR LEFT-OVER DATA AFTER ALL DEFINED
C  DATA HAS BEEN EXTRACTED.

C  NERR = RETURN CODE STEPED 1 FOR EACH ERROR.
C         (INPUT AND OUTPUT, INTEGER)

C*********************************************************************
      USE scrat1mod
      PARAMETER (MAXCHR=10000000)

      LOGICAL LQ
      CHARACTER*1 BLANK,COLON,QUOTE

      INCLUDE 'control.h'
!      INCLUDE 'scrat1.h'

      DATA BLANK/' '/,COLON/':'/,QUOTE/'"'/

      IF (MYPRC.GT.0) GO TO 5

      I=0
    1 I=I+1
      IF (I.GE.LAST) GO TO 5
      IF (A(I).NE.BLANK.AND.A(I).NE.COLON) THEN
         NERR=NERR+1
         LEVERR=2
         LQ=.TRUE.
         L2=I
         DO 2 J=I,LAST
         IF (A(J).EQ.QUOTE) LQ=.NOT.LQ
         IF ((A(J).EQ.COLON).AND.LQ) GO TO 3
    2    L2=J
    3    J2=L2
         IF (J2.GT.I+70) J2=I+70
         LEVERR=2
         IF (LEVELC) WRITE (NFOUT,4) (A(J),J=I,J2)
    4    FORMAT(/' ERROR # 116, UNDEFINED VARIABLE OR SYNTAX ERROR'/
     &      1X,75A1)
         I=L2
      ENDIF
      GO TO 1

    5 CONTINUE
      CALL SPREAD(1,NERR)

      END
C*********************************************************************
      SUBROUTINE GETGRDA (VNAM,KARY,NUMRET,NERR)
C*********************************************************************
C  ROUTINE READS A GRID-ELEMENT ARRAY.  DATA IS READ FOR ONLY THE
C  ELEMENTS ASSIGNED TO A PROCESSOR

C  VNAM   = VARIABLE NAME AND OPTIONAL UNITS (INPUT, CHARACTER*60).
C           THE NAME MUST BE TERMINATED WITH A BLANK OR THE LEFT
C           BRACKET OF A UNITS SPECFICATION.  THE NAME CAN NOT INCLUDE
C           EMBEDDED BLANKS.  UNITS, IF ANY, MUST BE ENCLOSED IN
C           BRACKETS [] AND IMMEDIATELY FOLLOW THE NAME.  BLANKS MAY BE
C           INCLUDED BETWEEN THE BRACKETS.  EXAMPLES:
C           'NX '   'P[psi]'   'HC[Btu/lb F]'

C  KARY  = ARRAY KEY (INPUT, INTEGER)
C        = 1 ==> BLOCK CENTER ARRAY
C        = 2 ==> BLOCK CORNER ARRAY
C        = 3 ==> MFMFE CORNER ARRAY

C  NUMRET = NUMBER OF VALUES READ TO THE ARRAY (OUTPUT, INTEGER)

C  NERR   = ERROR NUMBER STEPPED BY 1 FOR EACH DATA ERROR INCOUNTERED
C           (INPUT AND OUTPUT, INTEGER)

C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'blkary.h'
      INCLUDE 'control.h'
      INCLUDE 'readdat.h'
      CHARACTER*1 VNAM(*)
      INTEGER KARY,NUMRET,NERR
      CHARACTER*50 VNAM2
      INTEGER KARR,I,J,KE
      LOGICAL ISBINARY
      EXTERNAL GRDAIN,GRDBIN4,GRDBIN8,GRDBINI
      CHARACTER*80 :: INFILE

      INTEGER LASTNBLK,BINU
      COMMON /GRDBINCOM/LASTNBLK,BINU

C  PUT ARRAY NAME AND DATA IN COMMON

      CALL ARYDAT(VNAM,NTYPGA,ND4GA,KARR,KE)
      IF (KE.GT.0) GO TO 13
      KNDARY=KARY
      J=0
      DO 1 I=1,60
      VNAMGA(I)=VNAM(I)
      IF (VNAMGA(I).EQ.']'.OR.(VNAMGA(I).EQ.' '.AND.J.EQ.0)) GO TO 2
      IF (VNAMGA(I).EQ.'[') J=I
    1 CONTINUE

C  GET DATA VIA A WORK ROUTINE TO PROVIDE MEMORY FOR THE DATA

    2 NUMRGA=0
      KERRGA=0

! bag8 - check for binary input
      IF (J.GT.0) I=J-1
      WRITE(VNAM2,'(<I>A)')VNAM(1:I)
      ISBINARY=.FALSE.
      DO J=1,BINARY_INPUTS
        IF (TRIM(VNAM2).EQ.TRIM(INCLUDE_BINARY(J))) THEN
          ISBINARY=.TRUE.
          EXIT
        ENDIF
      ENDDO

      CALL ALLBLOCKS()
      IF (ISBINARY) THEN
        I4UTIL=KARY
!        IF (MYPRC.EQ.0) WRITE(*,*)'GETGRDA: ',TRIM(VNAM2),' (binary)'
        INFILE=TRIM(VNAM2)//'.bin'
        BINU=10
        OPEN(unit=BINU,file=TRIM(INFILE),
     &      form='binary',
     &      access='sequential',
     &      status='old')
        IF (MYPRC.EQ.0) WRITE(*,*)'Opening binary file: ', TRIM(INFILE)
        LASTNBLK=0
        IF (NTYPGA.EQ.1) THEN
          CALL CALLWORK(GRDBIN4,[2,KARR,N_I4U])
        ELSEIF (NTYPGA.EQ.2) THEN
          CALL CALLWORK(GRDBIN8,[2,KARR,N_I4U])
        ELSEIF ((NTYPGA.EQ.3).OR.(NTYPGA.EQ.4)) THEN
          CALL CALLWORK(GRDBINI,[2,KARR,N_I4U])
        ENDIF
        CLOSE(BINU)
      ELSE
!        IF (MYPRC.EQ.0) WRITE(*,*)'GETGRDA: ',TRIM(VNAM2),' (ascii)'
        CALL CALLWORK(GRDAIN,[1,KARR])
      ENDIF
      NUMRET=NUMRGA
      NERR=NERR+KERRGA

C  EXITS

      RETURN
   13 NERR=NERR+1
      END

C*********************************************************************
      SUBROUTINE GRDAIN (IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     & KEYOUT,NBLKA,VAL)
C*********************************************************************

C  WORK ROUTINE FOR GRID-ELEMENT ARRAY INPUT

C  VAL()  = ARRAY RETURNED (OUTPUT).

C*********************************************************************
C      INCLUDE 'msjunk.h'

      PARAMETER (NVN=60+4)

      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      CHARACTER*1 VNAMC(NVN),DIGITS(10)
      CHARACTER*2 VTYP(6)
      INCLUDE 'readdat.h'

      DATA VTYP/'R4','R8','I2','I4','L2','L4'/,
     & DIGITS/'0','1','2','3','4','5','6','7','8','9'/

C  BUILD VARIABLE NAME FOR A SPECIFIC GRID BLOCK

      DO 1 I=1,60
      L=I
      IF (VNAMGA(I).EQ.'['.OR.VNAMGA(I).EQ.' ') GO TO 2
    1 VNAMC(L)=VNAMGA(I)

    2 J=L
      IF (NBLKA.GT.99) THEN
         VNAMC(L)=DIGITS(1+NBLKA/100)
         L=L+1
      ENDIF
      IF (NBLKA.GT.9) THEN
         VNAMC(L)=DIGITS(1+MOD(NBLKA,100)/10)
         L=L+1
      ENDIF
      VNAMC(L)=DIGITS(1+MOD(NBLKA,10))
      L=L+1

      IF (VNAMGA(J).EQ.'[') THEN
         DO 3 I=J,60
         VNAMC(L)=VNAMGA(I)
         L=L+1
         IF (VNAMGA(I).EQ.']') GO TO 4
    3    CONTINUE
      ELSE
         IF (L.LE.NVN) VNAMC(L)=' '
      ENDIF

C  GET GRID BLOCK GLOBAL DIMENSIONS

    4 CALL BLKDIM(NBLKA,IDIMG,JDIMG,KDIMG,NERR)
      IF (NERR.NE.0) THEN
         KERRGA=KERRGA+1
         RETURN
      ENDIF
      IF (KNDARY.EQ.2.OR.KNDARY.EQ.3) THEN
         IDIMG=IDIMG+1
         JDIMG=JDIMG+1
         KDIMG=KDIMG+1
      ENDIF

C  PUT LOCAL-GLOBAL OFFSETS AND LOCAL DIMENSIONS IN /READAT/

      CALL BLKOFF(NBLKA,IGLT,JGLT,KGLT,NERR)
      IF (KNDARY.EQ.3) THEN
        IDIML=IDIM+1
        JDIML=JDIM+1
        KDIML=KDIM+1
      ELSE
        IDIML=IDIM
        JDIML=JDIM
        KDIML=KDIM
      ENDIF

C  GET DATA FOR THE GRID BLOCK AND RETURN

      NERR=0
      NBLKG=NBLKA
      IF (KNDARY.EQ.3) THEN
C         CALL GETVAL2(VNAMC,VAL,VTYP(NTYPGA),IDIMG,JDIMG,KDIMG,
C     & ND4GA,N,NERR)
      ELSE
      CALL GETVAL (VNAMC,VAL,VTYP(NTYPGA),IDIMG,JDIMG,KDIMG,ND4GA,N,
     & NERR)
      ENDIF
      NBLKG=0
      NUMRGA=NUMRGA+N
      IF (NERR.NE.0) KERRGA=KERRGA+1
      END

! bag8 - binary input routines
C*********************************************************************
      SUBROUTINE GRDBIN4(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     & KEYOUT,NBLK,VAL,KARY)
C*********************************************************************
C  bag8 - WORK ROUTINE FOR BINARY GRID-ELEMENT ARRAY INPUT
C       VAL()  = REAL*4 ARRAY RETURNED (OUTPUT).
C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'blkary.h'
      INCLUDE 'control.h'
      INCLUDE 'readdat.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,KARY
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*4 VAL(*)
      INTEGER I,J,K,L,IOFF,JOFF,KOFF,NERR,NX,NY,NZ,N,IDX
      INTEGER ID1,JD1,KD1,MP
      REAL*8, ALLOCATABLE :: BUF(:,:,:,:)
      INTEGER LASTNBLK,BINU
      COMMON /GRDBINCOM/LASTNBLK,BINU

! Skip ahead in the binary file if fault blocks do not belong
! to this processor.
      IF (LASTNBLK.NE.NBLK-1) THEN
      DO I=LASTNBLK+1,NBLK-1
        IF (KARY.EQ.1) THEN    ! cell-centered data
          NX=NXDIM(NBLK)
          NY=NYDIM(NBLK)
          NZ=NZDIM(NBLK)
        ELSE                   ! corner point data
          NX=NXDIM(NBLK)+1
          NY=NYDIM(NBLK)+1
          NZ=NZDIM(NBLK)+1
        ENDIF
        N=NX*NY*NZ*ND4GA
        ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
        READ(BINU) BUF
        DEALLOCATE(BUF)
      ENDDO
      ENDIF

! Read in current fault block to buffer
      IF (KARY.EQ.1) THEN
        NX=NXDIM(NBLK)
        NY=NYDIM(NBLK)
        NZ=NZDIM(NBLK)
        ID1=IDIM
        JD1=JDIM
        KD1=KDIM
        MP=0
      ELSE
        NX=NXDIM(NBLK)+1
        NY=NYDIM(NBLK)+1
        NZ=NZDIM(NBLK)+1
        ID1=IDIM+1
        JD1=JDIM+1
        KD1=KDIM+1
        MP=1
      ENDIF
      N=NX*NY*NZ*ND4GA
      ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
      READ(BINU) BUF

! Transfer data on this process to grid element array
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,NERR)
      DO L = 1,ND4GA
!      DO K = KL1,KL2+MP
!      DO J = JL1V(K),JL2V(K)+MP
!      DO I = IL1,IL2+MP
      DO K = 1,KD1
      IF ((K+KOFF.LT.1).OR.(K+KOFF.GT.NZ)) CYCLE
      DO J = 1,JD1
      IF ((J+JOFF.LT.1).OR.(J+JOFF.GT.NY)) CYCLE
      DO I = 1,ID1
      IF ((I+IOFF.LT.1).OR.(I+IOFF.GT.NX)) CYCLE
        IDX = I + (J-1)*ID1 + (K-1)*ID1*JD1 +(L-1)*ID1*JD1*KD1
        VAL(IDX)=BUF(I+IOFF,J+JOFF,K+KOFF,L)

! bag8 debug
!        WRITE(*,*)'NBLK,I,J,K,VAL=',NBLK,I+IOFF,J+JOFF,K+KOFF,VAL(IDX)
!        IF (VAL(IDX).EQ.0.D0) THEN
!          WRITE(*,*)'Zero value'
!          PAUSE
!        ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DEALLOCATE(BUF)

      NUMRGA=NUMRGA+N
      LASTNBLK=NBLK

      END

C*********************************************************************
      SUBROUTINE GRDBIN8(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     & KEYOUT,NBLK,VAL,KARY)
C*********************************************************************
C  bag8 - WORK ROUTINE FOR BINARY GRID-ELEMENT ARRAY INPUT
C       VAL()  = REAL*8 ARRAY RETURNED (OUTPUT).
C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'blkary.h'
      INCLUDE 'control.h'
      INCLUDE 'readdat.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,KARY
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      REAL*8 VAL(*)
      INTEGER I,J,K,L,IOFF,JOFF,KOFF,NERR,NX,NY,NZ,N,IDX
      INTEGER ID1,JD1,KD1,MP
      REAL*8, ALLOCATABLE :: BUF(:,:,:,:)
      INTEGER LASTNBLK,BINU
      COMMON /GRDBINCOM/LASTNBLK,BINU

! Skip ahead in the binary file if fault blocks do not belong
! to this processor.
      IF (LASTNBLK.NE.NBLK-1) THEN
      DO I=LASTNBLK+1,NBLK-1
        IF (KARY.EQ.1) THEN    ! cell-centered data
          NX=NXDIM(NBLK)
          NY=NYDIM(NBLK)
          NZ=NZDIM(NBLK)
        ELSE                   ! corner point data
          NX=NXDIM(NBLK)+1
          NY=NYDIM(NBLK)+1
          NZ=NZDIM(NBLK)+1
        ENDIF
        N=NX*NY*NZ*ND4GA
        ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
        READ(BINU) BUF
        DEALLOCATE(BUF)
      ENDDO
      ENDIF

! Read in current fault block to buffer
      IF (KARY.EQ.1) THEN
        NX=NXDIM(NBLK)
        NY=NYDIM(NBLK)
        NZ=NZDIM(NBLK)
        ID1=IDIM
        JD1=JDIM
        KD1=KDIM
        MP=0
      ELSE
        NX=NXDIM(NBLK)+1
        NY=NYDIM(NBLK)+1
        NZ=NZDIM(NBLK)+1
        ID1=IDIM+1
        JD1=JDIM+1
        KD1=KDIM+1
        MP=1
      ENDIF
      N=NX*NY*NZ*ND4GA
      ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
      READ(BINU) BUF

! Transfer data on this process to grid element array
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,NERR)
      DO L = 1,ND4GA
      DO K = 1,KD1
      IF ((K+KOFF.LT.1).OR.(K+KOFF.GT.NZ)) CYCLE
      DO J = 1,JD1
      IF ((J+JOFF.LT.1).OR.(J+JOFF.GT.NY)) CYCLE
      DO I = 1,ID1
      IF ((I+IOFF.LT.1).OR.(I+IOFF.GT.NX)) CYCLE
        IDX = I + (J-1)*ID1 + (K-1)*ID1*JD1 +(L-1)*ID1*JD1*KD1
        VAL(IDX)=BUF(I+IOFF,J+JOFF,K+KOFF,L)
        WRITE(*,*)'NBLK,I,J,K,VAL=',NBLK,I+IOFF,J+JOFF,K+KOFF,VAL(IDX)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DEALLOCATE(BUF)

      NUMRGA=NUMRGA+N
      LASTNBLK=NBLK

      END

C*********************************************************************
      SUBROUTINE GRDBINI(IDIM,JDIM,KDIM,LDIM,IL1,IL2,JL1V,JL2V,KL1,KL2,
     & KEYOUT,NBLK,VAL,KARY)
C*********************************************************************
C  bag8 - WORK ROUTINE FOR BINARY GRID-ELEMENT ARRAY INPUT
C       VAL()  = INTEGER ARRAY RETURNED (OUTPUT).
C*********************************************************************
      IMPLICIT NONE
C      INCLUDE 'msjunk.h'
      INCLUDE 'blkary.h'
      INCLUDE 'control.h'
      INCLUDE 'readdat.h'
      INTEGER IDIM,JDIM,KDIM,LDIM,IL1,IL2,KL1,KL2,NBLK,KARY
      INTEGER JL1V(KDIM),JL2V(KDIM),KEYOUT(IDIM,JDIM,KDIM)
      INTEGER VAL(*)
      INTEGER I,J,K,L,IOFF,JOFF,KOFF,NERR,NX,NY,NZ,N,IDX
      INTEGER ID1,JD1,KD1,MP
      REAL*8, ALLOCATABLE :: BUF(:,:,:,:)
      INTEGER LASTNBLK,BINU
      COMMON /GRDBINCOM/LASTNBLK,BINU

! Skip ahead in the binary file if fault blocks do not belong
! to this processor.
      IF (LASTNBLK.NE.NBLK-1) THEN
      DO I=LASTNBLK+1,NBLK-1
        IF (KARY.LE.1) THEN    ! cell-centered data
          NX=NXDIM(NBLK)
          NY=NYDIM(NBLK)
          NZ=NZDIM(NBLK)
        ELSE                   ! corner point data
          NX=NXDIM(NBLK)+1
          NY=NYDIM(NBLK)+1
          NZ=NZDIM(NBLK)+1
        ENDIF
        N=NX*NY*NZ*ND4GA
        ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
        READ(BINU) BUF
        DEALLOCATE(BUF)
      ENDDO
      ENDIF

! Read in current fault block to buffer
      IF (KARY.LE.1) THEN
        NX=NXDIM(NBLK)
        NY=NYDIM(NBLK)
        NZ=NZDIM(NBLK)
        ID1=IDIM
        JD1=JDIM
        KD1=KDIM
        MP=0
      ELSE
        NX=NXDIM(NBLK)+1
        NY=NYDIM(NBLK)+1
        NZ=NZDIM(NBLK)+1
        ID1=IDIM+1
        JD1=JDIM+1
        KD1=KDIM+1
        MP=1
      ENDIF
      N=NX*NY*NZ*ND4GA
      ALLOCATE(BUF(NX,NY,NZ,ND4GA),STAT=NERR)
      READ(BINU) BUF

! Transfer data on this process to grid element array
      CALL BLKOFF(NBLK,IOFF,JOFF,KOFF,NERR)
      DO L = 1,ND4GA
      DO K = 1,KD1
      IF ((K+KOFF.LT.1).OR.(K+KOFF.GT.NZ)) CYCLE
      DO J = 1,JD1
      IF ((J+JOFF.LT.1).OR.(J+JOFF.GT.NY)) CYCLE
      DO I = 1,ID1
      IF ((I+IOFF.LT.1).OR.(I+IOFF.GT.NX)) CYCLE
        IDX = I + (J-1)*ID1 + (K-1)*ID1*JD1 +(L-1)*ID1*JD1*KD1
        VAL(IDX)=BUF(I+IOFF,J+JOFF,K+KOFF,L)
      !  WRITE(*,*)'NBLK,I,J,K,VAL=',NBLK,I+IOFF,J+JOFF,K+KOFF,VAL(IDX)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DEALLOCATE(BUF)

      NUMRGA=NUMRGA+N
      LASTNBLK=NBLK

      END

C*********************************************************************
      SUBROUTINE DEFAULT (DEFUNT)
C*********************************************************************
C  ROUTINE SETS DEFAULT UNITS FOR THE NEXT (ONLY) CALL TO GETVAL OR
C  GETVALS.  ANY SPECIFICATION OF UNITS BY THE USER IN THE INPUT DATA
C  REPLACES THIS DEFAULT.

C  DEFUNT  = DEFAULT EXTERNAL UNITS (INPUT, CHARACTER*1)
C            TERMINATE WITH ]

C*********************************************************************
      CHARACTER*1 DEFUNT(*)
      CHARACTER*60 UNTD

      INCLUDE 'readdat.h'

      EQUIVALENCE (UNTD,UNTDEF(1))

      ISUNTD=.TRUE.
      UNTD=' '
      DO 1 I=1,60
      UNTDEF(I)=DEFUNT(I)
      IF (DEFUNT(I).EQ.']') RETURN
    1 CONTINUE
      END
C*********************************************************************
      SUBROUTINE GETBLK (VNAM,BUF,MAXBUF,NDIM,LOC,LEN,NERR)
C*********************************************************************
C  READS A SUBSCRIPTED BLOCK TEXT VARIABLE

C  VNAM   = VARIABLE NAME (INPUT, CHARACTER*60)
C           THE NAME MUST BE TERMINATED WITH A BLANK

C  BUF()  = SPACE FOR TEXT RETURNED (OUTPUT, CHARACTER*1)

C  MAXBUF = MAX CHARACTERS IN BUF() (INPUT, INTEGER)

C  NDIM   = DIMENSION OF THE BLOCK TEXT VARIABLE (INPUT, INTEGER)
C           SET NDIM = 0 IF THE VARIABLE IS NOT DIMENSIONED

C  LOC()  = LOCATIONS OF THE BLOCKS IN BUF() (OUTPUT, INTEGER)
C           DIMENSION NDIM

C  LEN()  = LENGTHS OF THE BLOCKS IN BUF() (OUTPUT, INTEGER)
C           DIMENSION NDIM

C  NERR   = ERROR NUMBER STEPPED BY 1 IF A DATA ERROR IS INCOUNTERED
C           (INPUT AND OUTPUT, INTEGER)

C  NOTE:  IF BLOCK I IS NOT READ THEN LOC(I) AND LEN(I) ARE NOT CHANGED
**********************************************************************

      INTEGER LOC(*),LEN(*)
      CHARACTER*1 VNAM(*),BUF(*)

      INCLUDE 'readdat.h'

      N=NDIM
      IF (N.LT.1) N=1
      DO 1 I=1,N
    1 LENBLK(I)=0
      NERRO=NERR
      CALL GETVALS(VNAM,BUF,'BT',NDIM,0,0,MAXBUF,NN,NERR)
      IF (NERR.NE.NERRO) RETURN
      DO 2 I=1,N
      IF (LENBLK(I).GT.0) THEN
         LOC(I)=LOCBLK(I)
         LEN(I)=LENBLK(I)
      ENDIF
    2 CONTINUE
      END
C*********************************************************************
      SUBROUTINE SETNBG ()
C*********************************************************************

C  SILLY LITTLE ROUTINE TO INITIALIZE NBLKG

C*********************************************************************
      INCLUDE 'readdat.h'
      NBLKG=0
      NOERASE=.FALSE.
      END
C*********************************************************************
      SUBROUTINE TBLKOUT (BLK,LB)
C*********************************************************************

C  PRINTS A BLOCK TEXT STRING

C  BLK() = BLOCK TEXT STRING (INPUT, CHARACTER*1)

C  LB = LENGTH OF BLK (INPUT, INTEGER)

C*********************************************************************
      INCLUDE 'control.h'
      CHARACTER*1 BLK(*),RECEND

      RECEND=CHAR(30)
      I1=1
      DO 1 I=1,LB
      IF (BLK(I).EQ.RECEND.OR.(I.EQ.LB.AND.I1.LT.LB)) THEN
         WRITE(NFOUT,2) (BLK(J),J=I1,I-1)
         I1=I+1
      ENDIF
    1 CONTINUE
    2 FORMAT(80A1)

      END
