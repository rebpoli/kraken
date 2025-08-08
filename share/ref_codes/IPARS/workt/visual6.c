// -----------------------------------------------------------------
// file: visual6.c
//
// Corner point visualization output in the Tecplot format for
// IPARS framework
//
// Corner point variables so that interpolation.
//-----------------------------------------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "cfsimple.h"
#include "visualc.h"

// ----------------------------------------------------------------
// the actual unstructured grid routine
//

void _vis_corner_point
(
 const _F_INTEGER * const IDIM,
  const _F_INTEGER * const JDIM,
  const _F_INTEGER * const KDIM,
  const _F_INTEGER * const LDIM,
  const _F_INTEGER * const IL1,
  const _F_INTEGER * const IL2,
  const _F_INTEGER * const JLV1,
  const _F_INTEGER * const JLV2,
  const _F_INTEGER * const KL1,
  const _F_INTEGER * const KL2,
  const _F_INTEGER * const KEYOUT,
  const _F_INTEGER * const NBLK,
 // fixed parameters
 const int flag,
 const _F_REAL_8 * const xc,
 const _F_REAL_8 * const yc,
 const _F_REAL_8 * const zc,
 const int nscl, const int nvec,
 _F_REAL_8 * * r8_scl_list,
 _F_REAL_8 * * r8_vec_list
 )
{
  const int narg=nscl+nvec;
  const int maxblk = (*IDIM)*(*JDIM)*(*KDIM);
  const int maxpls = (*IDIM +1)*(*JDIM +1)*(*KDIM +1);
  const int   kinc = (*IDIM) * (*JDIM) ;
  const int   jinc = (*IDIM) ;
  const int   kbeg = (*KL1) ;
  const int   kend = (*KL2) ;
  const int   ibeg =  *IL1 ;
  const int   iend =  *IL2 ;
  const int   nblk = *NBLK-1;
  //const int _Myblk = nblk;

  //
  FILE *fp; // general file handle, first used for INITfile, next for ZONEfile
  int gi;

  if(flag != 6) return ;

  // DEBUG
  // printf("\n Unstruct routine with visflag=%d narg=%d\n",flag,narg);

  if (_Vis_cnt[nblk]++ ==0) {

    char INITfilename[50];
    _cr_INIT_fname(INITfilename,flag,nblk);

    fp=fopen(INITfilename,"w");
    if (fp != NULL) {
      int k,j,i,iarg,iblk=0,ic=0;
      _F_INTEGER igoff,jgoff,kgoff,nofferr;
      _F_INTEGER idgoff,jdgoff,kdgoff;

      _INFO_print(0,INITfilename,nblk);

      //
      // first count the gridblocks : necessary to determine
      // everything else
      //

      for(k=kbeg;k<=kend;k++) {
	const int jbeg =  ( JLV1[k-1] );
	const int jend =  ( JLV2[k-1] );	   	    	
	for ( j = jbeg ; j <= jend ;j++ ) {
	  for ( i= ibeg; i <= iend ; i++) {	      	
#define gindex(i,j,k)   ( ( (k)-1 )*kinc+ ( (j)-1)*jinc + ( (i) -1) )
	    const int index = gindex(i,j,k);	
#undef gindex	
	    if ( KEYOUT[index] == 1 ) 	
	      iblk++;	
	  } // end i
	} // end j
      } // end k
      _Ngblk[nblk]=iblk;

      //
      // print the header to the file
      //
      {
	int iarg;
	
	fprintf(fp,"title=init.%d\n",flag);
	fprintf(fp,"variables=x,y,z");

	// ignore the vector variables

	for(iarg=0;iarg<nscl;iarg++) {
	  if(_Vis_pars[iarg].name!=NULL)
	    fprintf(fp,", %s",_Vis_pars[iarg].name);
	  else
	    fprintf(fp,", v%d",iarg+1);
	}
      }

      //
      // print the xyz information
      //

      myblkoff(NBLK,&igoff,&jgoff,&kgoff,&nofferr);
      idgoff=igoff;jdgoff=jgoff;kdgoff=kgoff;

      // modify offsets to include offset of (0,0,0) block
      igoff+=_Mxrecxp*nblk;jgoff+=_Myrecyp*nblk;kgoff+=_Mzreczp*nblk;
      idgoff+=(_Mxrecxp-1)*nblk;jdgoff+=(_Myrecyp-1)*nblk;
      kdgoff+=(_Mzreczp-1)*nblk;

      // DEBUG
      /** printf("\n Offsets are :%d %d %d err=%d\n",
	     igoff,jgoff,kgoff,nofferr);
	     **/
      fprintf(fp,
	      "\nzone t=init, n=%d, e=%d, et=brick F=fepoint\n",
	      8*_Ngblk[nblk],_Ngblk[nblk]);
    	
      for(k=kbeg;k<=kend;k++) {
	const int jbeg =  ( JLV1[k-1] );
	const int jend =  ( JLV2[k-1] );	   	    	
	for ( j = jbeg ; j <= jend ;j++ ) {
	  for ( i= ibeg; i <= iend ; i++) {	      	
	    //
#define gindex(i,j,k)   ( ( (k)-1 )*kinc+ ( (j)-1)*jinc + ( (i) -1) )
	    // could be done more efficiently but this is readable
	
	    const int index = gindex(i,j,k);
	
	    if ( KEYOUT[index] == 1 ) {	

	      // printf x,y,z coordinates of the 8 nodes and
	      // print dummy=0.0 values of the vis variables

	      gi = index;
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);	     	
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i+1,j,k);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i,j+1,k);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i+1,j+1,k);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);	     	
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i,j,k+1);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i+1,j,k+1);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);	     	
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i,j+1,k+1);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);	     	
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	      gi = gindex(i+1,j+1,k+1);
	      fprintf(fp,"%g %g %g\n",xc[gi],yc[gi],zc[gi]);	     	
	      for(iarg=0;iarg<narg;iarg++)fprintf(fp,"%g ",0.0);
	      fprintf(fp,"\n");

	    } // end if KEYOUT
#undef gindex	
	  } // end i
	} // end j
      } // end k

      //
      // write the connectivity data
      //

      ic=0;
      for(k=kbeg;k<=kend;k++) {
	const int jbeg =  ( JLV1[k-1] );
	const int jend =  ( JLV2[k-1] );	   	    	
	for ( j = jbeg ; j <= jend ;j++ ) {
	  for ( i= ibeg; i <= iend ; i++) {	      	
#define gindex(i,j,k)   ( ( (k)-1 )*kinc+ ( (j)-1)*jinc + ( (i) -1) )
	    const int index = gindex(i,j,k);	
	
	    if ( KEYOUT[index] == 1 ) {
	      fprintf(fp,"%d %d %d %d %d %d %d %d\n",
		      ic+1,ic+2,ic+4,ic+3,ic+5,ic+6,ic+8,ic+7);
	      ic+=8;
	    } // end if keyout
#undef gindex
	  } // end i
	} // end j
      } // end k

      fclose(fp);

      // DEBUG
      //printf("\nVIS_CONST_UNSTRUCT info: %d gridblocks\n",
      //iblk);

    } // if file existed

  } // end if first time : Vis_cnt=0


  // ----------------------------------------------
  // print the values of variables to the ZONE file
  //
  {
    char ZONEfilename[50];
    _cr_ZONE_fname(ZONEfilename,flag,nblk);

    fp=fopen(ZONEfilename,"w");
    if(fp !=NULL) _INFO_print(1,ZONEfilename,nblk);
    else
      {
	printf("\n Problems with file: %s\n",ZONEfilename);exit(-1);
      }
  }

  if (fp != NULL) {
    int lvar,i,j,k;

    static char ZONEname [50];
    sprintf(ZONEname,"zone%d.%d.%d",nblk,_Myprc,_Nstep);

    // give info about the zone data

    fprintf(fp,
     "\nzone t=%s, n=%d, e=%d, et=brick, F=feblock,d=(1,2,3,FECONNECT)\n",
	    ZONEname,8*_Ngblk[nblk],_Ngblk[nblk]);

    // print the values of the scalar variables as cell-centered data
    // 8 equal values for a block will create a piecewise constant field

    for(lvar=0;lvar<nscl;lvar++) {
      const int ldimoff=(_Vis_pars[lvar].ldim-1)*maxblk;

      // printf("\nprinting for scalar variable %d\n",lvar);

      fprintf(fp,"\n\n");

      for(k=kbeg;k<=kend;k++) {
	const int jbeg =  ( JLV1[k-1] );
	const int jend =  ( JLV2[k-1] );	   	    	
	for ( j = jbeg ; j <= jend ;j++ ) {
	  for ( i= ibeg; i <= iend ; i++) {	      	
#define gindex(i,j,k)   ( ( (k)-1 )*kinc+ ( (j)-1)*jinc + ( (i) -1) )
	    const int index = gindex(i,j,k);	
#undef gindex	    	
	    if ( KEYOUT[index] == 1 ) {	
	      _F_REAL_8 *ptr=& ( r8_scl_list[lvar][ldimoff] );
	      _F_REAL_8 scalval=ptr[index];	
	      int repeat;

	      for(repeat=0;repeat<8;repeat++)
		fprintf(fp,"%g ",scalval);
	      fprintf(fp,"\n");

	    } // if keyout
	  } // if i
	} // if j
      } // if k
    } // lvar
	
    // NOTE : the vector variables need to be implemented!

    fprintf(fp,"\n");
    fclose(fp);

  } // ZONEfile open

}

