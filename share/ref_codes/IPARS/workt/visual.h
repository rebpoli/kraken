c ------------------------------------------------------------------
C  visual.h - visualization parameters and DATA for IPARS framework
c MPeszynska, 5/27/98	initial version
c MPeszynska, 9/25/98   changed some parameters to be keyword dependent
c                       added provisions for Multi Model
c Ben Ganis, Gergina Pencheva 2/7/13   added VIS_DIR
c Gergina Pencheva 10/12/15  added flag for ascii/binary vtk output
c -----------------------------------------------------------------
      INTEGER MAXVISVARS, MAXMODELS, MAXVISFLAG
      PARAMETER (MAXVISVARS = 30, MAXMODELS =19)
      PARAMETER (MAXVISFLAG = 11)

      CHARACTER*10
     &     VIS_VARNAMES (MAXVISVARS)
      CHARACTER*10  VIS_IPARS_NAMES (MAXVISVARS)
      CHARACTER*50 VIS_FNAME
      CHARACTER*10  VIS_FSTYLE
      CHARACTER*50 VIS_DIR

C      LOGICAL  IFPV3
C      CHARACTER*1 PV3_KEY(MAXVISVARS)
C      REAL*4      PV3_LIM1(MAXVISVARS),PV3_LIM2(MAXVISVARS)
C      COMMON /VISPV3/ PV3_KEY,PV3_LIM1,PV3_LIM2,IFPV3

      LOGICAL   IFPERM,IFWELL,IFVISINIT
      LOGICAL   VIS_BINARY,TECEXISTS,ALLZONES
      INTEGER	VISFLAG,N_VISFLAG,N_VIS_BINARY,
     &    N_VISDUM,VIS_NVARS,N_VIS_NVARS,N_LFACES,
C --- SAUMIK,BGANIS
     &    VIS_SCL_POROHEX,VIS_SCL_FLOW(MAXMODELS),STARTMODACT,
     &    N_VIS_SCL_POROHEX,N_VIS_SCL_FLOW(MAXMODELS),
C --- SAUMIK,BGANIS
     &    VIS_SCL,VIS_VEC,N_VIS_SCL,N_VIS_VEC,VIS_NVEC,
     &    N_VIS_XREC,N_VIS_YREC,N_VIS_ZREC,
     &    N_VIS_dXREC,N_VIS_dYREC,N_VIS_dZREC,
     &    VIS_NAMELEN, N_VIS_VARNAMES,
     &    N_VIS_VARS(MAXVISVARS,MAXMODELS), VIS_OFFSETS(MAXVISVARS),
     &    VIS_VAL_NODAL(MAXVISVARS)

      COMMON /visualc/ VIS_VAL_NODAL, VISFLAG, N_VISFLAG,
     &    N_VISDUM,IFPERM,IFWELL,IFVISINIT,
     &    VIS_BINARY, TECEXISTS, ALLZONES, N_VIS_BINARY,
     &    VIS_NVARS,N_VIS_NVARS,N_LFACES,
C --- SAUMIK,BGANIS
     &    VIS_SCL_POROHEX,VIS_SCL_FLOW,STARTMODACT,
     &    N_VIS_SCL_POROHEX,N_VIS_SCL_FLOW,
C --- SAUMIK,BGANIS
     &    VIS_SCL,VIS_VEC,N_VIS_SCL,
     &    N_VIS_VEC,VIS_NVEC,N_VIS_XREC,N_VIS_YREC,N_VIS_ZREC,
     &    N_VIS_dXREC,N_VIS_dYREC,N_VIS_dZREC,VIS_NAMELEN,
     &    N_VIS_VARNAMES,N_VIS_VARS,VIS_OFFSETS,VIS_VARNAMES,
     &    VIS_FNAME,VIS_FSTYLE,VIS_IPARS_NAMES,VIS_DIR

      LOGICAL FIRSTVAR
      COMMON /VISFIRST/ FIRSTVAR

c VISFLAG - flags for types of visualization output
c
c       Old Tecplot formats (uses interpolation, requires postprocessing)
c
c          = 0     no visualization output
c          = 1     structured output, brick models
c                  (only works on one processor for one fault block)
c                        cell-centered pressures (shifted to
c                            corners) to tec.proc.tstep
c                        velocities (same as above, cell centered
c			     values are really averages of the
c			     adjacent faces) to vel.proc.tstep
c
c          = 2     unstructured output, brick models
c                        discontinuous pressure/sat data from cell-centers
c                        cell-centered values mapped to all the 8 vertices
c
c          = 3     unstructured output, brick models
c                        nodal values are interpolated from cell-centers
c                        nodal velocities: computed through interpolation
c
c          = 4     rectangle output, brick models
c
c          = 5     binary unstructured output, brick models
c                  (each variable in separate file)
c                        nodal values are interpolated from cell-centers
c                        nodal velocities: computed through interpolation
c
c          = 6     corner point output, hexahedral models
c                  (Not currently working for hexahedral MFMFE models)
c
c bag8 - VTK formats
c
c          = 7     bricks
c          = 8     hexahedra
c
c bag8 - Modern TecPlot formats
c          (cell-centers without interpolation, no postprocessing)
c
c          = 9     bricks or hexahedra
c          = 11    polyhedra
c
c bag8 - Both modern TecPlot and VTK
c
c          = 10    bricks or hexahedra (functions as both 7 & 9, or 8 & 9)

c VIS_BINARY (used when VISFLAG = 7,8,9)
c          = TRUE    output files are binary (default)
c          = FALSE   output files are ascii

c N_VISDUM : grid array number of dummy visualization grid array:
c this one holds values 0.0 used when wrong variable name is set in input
c file

c ifperm : flag which speicfies whether the initial outpout (at time zero)
c     should be created : the call is from inside of IVPARM, before TRANC1

c VIS_NVARS, N_VIS_NVARS: number of variables to be visualized - it really
c is the number of pointers passed and different names, for example
c if you use DUNK(1) and DUNK(27) then they appear twice
c and so the VIS_NVARS=2 not 1
c
c VIS_SCL: how many scalar variables
c VIS_NVEC: how many vector variables
c VIS_VEC: total of vector variables, equal VIS_NVEC*3 , 3 =XYZ components
c
c N_VIS_XREC,N_VIS_YREC,N_VIS_ZREC - IPARS pointers to
c           the xrec,yrec,zrec data
c
c VIS_VARNAMES = (character) array of variable names for output
c VIS_IPARS_NAMES = (character) array of IPARS names for variables
c
c VIS_NAMELEN = length of each of the characters in VIS_VARNAMES
c
c VIS_FNAME =
c 	the file name for the tecplot output, default
c 	is TEC.faultblock.processor.visflag.init
c 	or TEC.faultblock.processor.visflag.timestep
c zone name is constructed similarly
c VIS_FSTYLE =
c 	flag that tells the file format for output
c 	currently available PC, NO_PC (default)
c
c N_VIS_VARS = contains IPARS pointers to the variables for diff. models
c VIS_OFFSETS = contains offsets for the variables
c

C FIRSTVAR (logical) is .true. if this is the first time a translation
c	routine (TVIS_TRANSL etc.) is called, .false. otherwise.
c  	May be used for translation of some extra variables for models
c









