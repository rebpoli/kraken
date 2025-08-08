C  SBLKC.H - MULTIBLOCK DATA FOR DUAL APPROXIMATION INTERFACE

C  CODE HISTORY:

C  JOHN WHEELER      2/17/99    ALPHA CODE
C  JOHN WHEELER      3/14/99    MEMORY MANAGEMENT FOR INTERFACE BUFFERS
C  JOHN WHEELER      5/28/99    MULTIMODEL CAPABILITY
C  RICK DEAN         1/07/02    ADDED SPECIFICATIONS FOR ALL VARIABLES
C  SUNIL G. THOMAS   -/--/--    MODS FOR COUPLING FLOW TO TRCHEM AND
C                               VECTOR VIS MODS FOR EVMFEM
C  BEN GANIS         12/2/15    LOCAL FLUX IMPLEMENTATION FOR MULTIBLOCK
C                               WITH GENERAL HEXAHEDRA

C*********************************************************************


      REAL*8 TFINS, FACECP, AREAI, AREAT, UNKMAP
      REAL*8 TDDFINS0,TDDFINS
      REAL*4 COFINF
      INTEGER NFACES, NIBLKS, NAEB, IIEBS, NFIEBS, IJKS, NCGES, ICGES,
     & NFICGES, JBLOCK, NKDIRS, KDIRS, LIBUF, NPAI, NERRI, NPSRI, NEIFS,
     & NBSRI, NESNDI, NFESR, NUMBUFU, N_BUFIF, N_BUFDIM, IESNDI,
     & KFESR, NIEBS, NERECI, IERECI
      INTEGER TRDDIJKT
      INTEGER IJKT(3,90000)

      COMMON /SBLKC/ TFINS(90000), AREAI(90000), FACECP(6,11),
     & AREAT(11),UNKMAP(3,3,19,19),
     & COFINF(90000,3,3),
     & NFACES, NIBLKS(2,11), NAEB(90000), IIEBS(10),
     & NFIEBS, IJKS(3,40000), NCGES(40000), ICGES(40000), NFICGES,
     & JBLOCK(90000), NKDIRS, KDIRS(90000), LIBUF(90000),
     & NPAI(10), NERRI, NPSRI(10,10), NEIFS(11),
     & NBSRI(10,10), NESNDI(10,10), NFESR, NUMBUFU,
     & N_BUFIF, N_BUFDIM, IESNDI(10,10), KFESR(90000),
     & NIEBS(10), NERECI(10,10), IERECI(10,10),
     & IJKT

      COMMON /STRDDBLKC/ TRDDIJKT(3,90000), TDDFINS0(90000),
     &                   TDDFINS(90000)


! bag8 - local flux
      INTEGER EVFEM_MAP
      COMMON /EVFEM_MAP_COM/EVFEM_MAP

C*********************************************************************

C From dual size file:

C  MXGEI  = Max total number of interface grid elements per processor
C  MXFELE = Max total number of interface couplings per processor

C From frame size file:

C  MXFACE = Max total number of interfaces
C  MXMPP  = Max number of adjacent processors including self
C  MXBLKS = Max number of fault blocks
C  MXNUMEQ = Max number of equations per element
C  MXMOD = Max number of physical models in framework

C*********************************************************************

C  NFACES    = Number of interfaces total.

C  NIBLKS(,n) = A and B fault block numbers for interface n

C  FACECP(,n) = X,Y,Z data in the coordinate systems of the two fault
C               blocks for a common point.  Order: Xa Ya Za Xb Yb Zb

C  NEIFS(n) = Number of grid element interactions on interface n

C  NIEBS(m) = Number of grid elements on any interface for fault block m
C             (current processor only).

C  IIEBS(m) = Packing offset of grid elements on any interface for fault
C             block m (current processor only) (associated with k index
C             in definitions below).

C  NFIEBS = Next free k index (initialize to 1).

C  IJKS(3,k) = Local I,J,K of grid element k on an interface (source block)

C  NCGES(k) = Number of grid elements connected to grid element k through
C             the interface

C  ICGES(k) = Packing offset of grid elements connected to grid element k
C             through the interface (associated with j index in
C             definitions below)

C  NFICGES = Next free j index (initialize to 1).

C  JBLOCK(j) = Fault block number of element j (target block)

C  IJKT(3,j) = Global I,J,K of grid element j on an interface (target block)
C              Used only during initialization.

C  TRDDIJKT(3,j) = Global I,J,K of grid element j on an interface (target block)
C              Reserved for diffusion-dispersion problem.

C  AREAI(j)  = Area of interfacial element intersection.

C  TFINS(j)  = 2 * A    / (D   / K   + D   / K  )
C                   AkBj    Ak    Ak    Bj    Bj

C  COFINF(j,M,N) = Interaction coefficient for equation m and variable n
C                = Derivative of Q      WRT V
C                                 MAkBj      NBj

C  NKDIRS = Number of active diretion keys in KDIRS().  Normally six.

C  KDIRS(j) = Direction key for element j relative to element k

C           =  1 ==> X + 1

C           =  2 ==> Y + 1

C           =  3 ==> Z + 1

C           =  4 ==> X - 1

C           =  5 ==> Y - 1

C           =  6 ==> Z - 1

C           =  Additional offsets for higher order approximations.

C  LIBUF(j) = Receive buffer index of data from across the interface
C             (ie. first index in BUFIF4(,) or BUFIF8(,) )

C  NPAI(m) =  Number of processors that interact with the current processor
C             and fault-block m cross interfaces.  If the current processor
C             has grid elements on both sides of an interface, include it in
C             the count.

C  NPSRI(i,m) = Processor to send/receive message i for fault-block m.

C  NBSRI(i,m) = Fault block to send/receive message i for fault-block m.

C  NESNDI(i,m) = Number of grid elements in send message i for fault-block n.

C  IESNDI(i,m) = Packing offset of grid elements in send message i for
C                fault-block m (associated with ii index in definition below).

C  NERECI(i,m) = Number of grid elements in receive message i for fault-block n.

C  IERECI(i,m) = Packing offset of grid elements in receive message i for
C                fault-block m.

C  NFESR = Next free buffer index (ii) (initialize to 1).

C  KFESR(ii) = k of grid element in a send message

C  NERRI = Interface error count.

C  N_BUFIF = Memory management index for BUFIF8 and BUFIF4 which are
C            available to physical models after initialization.
C            Dimensioning:  BUFIF8(NFESR,34) or BUFIF4(2*NFESR,34)

C  N_BUFDIM = Memory management index for first dimension of BUFIF8(), NFESR

C  NUMBUFU  = Utility variable for passing a buffer number into interface
C             work routines

C  NAEB(j)  = k of grid element.  Used by framework only during
C             initialization.

C  AREAT(n) = Total area of interface n.

C  UNKMAP(nv1,nv2,nm1,nm2) = Matrix that maps the primary variables in model
C                            nm1 to to the variables in model nm2.

C bag8

C  EVFEM_MAP = choice of hard coded global mapping when EVFEM_HEX=2
C    (See: subroutine MAP in infcomm.df)

