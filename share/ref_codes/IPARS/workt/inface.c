// inface.c - C routines for Dual Approximatin Interface

// HISTORY:

// JOHN WHEELER      12/24/98   ORIGINAL BETA CODE
// JOHN WHEELER       3/14/99   MEMORY MANAGEMENT FOR INTERFACE BUFFERS

#include "memory.h"

// ROUTINE DECLARATIONS

FORTSUB callface_(PINT4 nif,PINT4 na,PINT4 nb);
FORTSUB alcibuf_(PINT4 numibuf,PINT4 sizibuf,PINT4 arynum,
        PINT4 err);
FORTSUB realcibuf_(PINT4 numibuf,PINT4 sizibuf,PINT4 arynum,
        PINT4 err);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FORTSUB callface_ (PINT4 nif, PINT4 nai, PINT4 nbi)
{
// *******************************************************************

// Calls SETINF routine with memory data.  Used during setup of the
// interface system.

// nif = Interface number (input)

// nai, nbi = Fault block numbers (input)

// *******************************************************************
int na, nb;

extern FORTSUB setinf_(PINT4 ida,PINT4 jda,PINT4 kda,PINT4 il1a,
       PINT4 il2a,PINT4 jl1a,PINT4 jl2a,PINT4 kl1a,PINT4 kl2a,PINT4 keya,
       PINT4 idb,PINT4 jdb,PINT4 kdb,PINT4 il1b,
       PINT4 il2b,PINT4 jl1b,PINT4 jl2b,PINT4 kl1b,PINT4 kl2b,PINT4 keyb,
       PINT4 nblka,PINT4 nblkb,PINT4 nif);

na = (*nai) - 1;
nb = (*nbi) - 1;

if (myelem[na] <= 0) return;

setinf_ (&(idim[na]),&(jdim[na]),&(kdim[na]),&(iloc1[na]),
   &(iloc2[na]),jloc1[na],jloc2[na],&(kloc1[na]),&(kloc2[na]),keyout[na],
   &(idim[nb]),&(jdim[nb]),&(kdim[nb]),&(iloc1[nb]),
   &(iloc2[nb]),jloc1[nb],jloc2[nb],&(kloc1[nb]),&(kloc2[nb]),keyout[nb],
   nai,nbi,nif);

return;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FORTSUB alcibuf_ (PINT4 numibuf, PINT4 sizibuf, PINT4 arynum,
        PINT4 err)
{
// *******************************************************************

// Alocates buffer space for the dual approximation block interface
// option.  Called only by the IPARS framework.

// numibuf = Number of buffers to be allocated (input).
//           Same number for all processors.

// sizibuf = Number of REAL*8 locations in one buffer (input).
//           May vary from one processor to the next.

// arynum = Array number (output)

// err = Error number (output)
//     = 0 ==> no error
//     = 464 ==> too many arrays
//     = 462 ==> insufficient memory

// *******************************************************************
int na, nb;
//long nm; //sgt
int nm;
PADD varb;

for (na = 0; (na < MAXARY) && (arytyp[na] != 0); na++);
if (na == MAXARY)
   {
   *err = 464L;
   return;
   }

arytyp[na] = 1;
aryfrm[na] = 3;

nm = (*numibuf) * (*sizibuf) * sizeof(double);
if ((varb = (PADD) malloc(nm)) == NULL)
   {
   *err = 462L;
   return;
   }

for (nb = 0; nb < numblks; nb++) aryadd[nb][na] = varb;

if (numarys == na) numarys++;
*arynum = na;
*err = 0L;
return;
}

// bag8
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FORTSUB realcibuf_(PINT4 numibuf,PINT4 sizibuf,PINT4 arynum,
        PINT4 err)
{
// *******************************************************************

// Alocates buffer space for the dual approximation block interface
// option.  Called only by the IPARS framework.

// numibuf = Number of buffers to be allocated (input).
//           Same number for all processors.

// sizibuf = Number of REAL*8 locations in one buffer (input).
//           May vary from one processor to the next.

// arynum = Array number (output)

// err = Error number (output)
//     = 0 ==> no error
//     = 464 ==> too many arrays
//     = 462 ==> insufficient memory

// *******************************************************************
int na, nb;
//long nm; //sgt
int nm;
PADD varb;

na = *arynum;

free(aryadd[0][na]);
for (nb = 0; nb < numblks; nb++) {
  aryadd[nb][na] = NULL;
}

nm = (*numibuf) * (*sizibuf) * sizeof(double);
if ((varb = (PADD) malloc(nm)) == NULL)
   {
   *err = 462L;
   return;
   }

for (nb = 0; nb < numblks; nb++) aryadd[nb][na] = varb;

*err = 0L;
return;
}
