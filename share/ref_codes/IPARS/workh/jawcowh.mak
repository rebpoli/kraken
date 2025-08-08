# Makefile Executive for IPARS on longhorn (ethernet) with single phase
# and hydrology models, lsor, dual, and mpi

# Copy this file to ipars/work/makefile then edit it to uncomment the
# appropriate lines

# make                 Builds production program
# make clean           Deletes work files

########################## Unix/Dos Controls #########################

# Define the slash for file names
S=/

# Define the target file name
EXENAM=ipars

# Define the object file extension
O=.o

# Define the copy instruction
COPY=cp

############################## Misc ##################################

default:$(EXENAM)

SIZE=.$(S)size
SIZDAT=jawcowh.siz
SETSIZE=echo $(SIZDAT) $@ > ech
SETSIZE1=echo $(SIZDAT) $< $@ > ech
SIZEIT=$(SIZE) < ech

MAKDIR=..$(S)make$(S)

WORK=.

COPYIT=$(COPY) $? $(WORK)

.SUFFIXES:
.SUFFIXES: .o .f .c .C .h .df .dc .dh .dC .in .obj

####################### Framework Include Files ######################

include $(MAKDIR)parall.mak
include $(MAKDIR)frame.mak

###################### Linear Solver (one only) ######################

include $(MAKDIR)lsor.mak
# include $(MAKDIR)gmres.mak
# include $(MAKDIR)petsc.mak
# include $(MAKDIR)pcg.mak 

####################### Multi-block (one only) #######################

# include $(MAKDIR)mortar.mak           
include $(MAKDIR)dual.mak

############################ Graphics ################################

include $(MAKDIR)visual.mak           

################### Physical Model Include Files #####################

include $(MAKDIR)hydroi.mak           
# include $(MAKDIR)hydroe.mak           
# include $(MAKDIR)blacki.mak           
# include $(MAKDIR)singlee.mak
# include $(MAKDIR)singlei.mak           
# include $(MAKDIR)comp.mak             Out of date

##################### Combine Object/Lib Files #######################

MODELOBJ = $(HYDROIOBJ) $(HYDROEOBJ) $(SINGLEIOBJ) $(SINGLEEOBJ) \
  $(BLACKIOBJ) $(COMPOBJ)

OBJS = $(FRAMEOBJ) $(MODELOBJ) $(BLOCKSOBJ) $(SOLVEOBJ) $(GRAPHOBJ) $(PARALOBJ)

MODELLIB = $(HYDROILIB) $(HYDROELIB) $(SINGLEILIB) $(SINGLEELIB) \
  $(BLACKILIB) $(COMPLIB)

LIBS = $(FRAMELIB) $(BLOCKSLIB) $(MODELLIB) $(SOLVELIB) $(GRAPHLIB) $(PARALLIB)

########### Machine and Compiler Include File (one only) #############
# This include file must be the last one

# include $(MAKDIR)sp2.mak
# include $(MAKDIR)lnx.mak
# include $(MAKDIR)beo.mak
# include $(MAKDIR)x2.mak
include $(MAKDIR)x2eth.mak
# include $(MAKDIR)t3e.mak              Untested
# include $(MAKDIR)rs6k.mak             Untested
