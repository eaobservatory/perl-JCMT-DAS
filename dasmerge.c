#ifdef DEBUG
#include <stdio.h>
#endif

#include "dasmerge.h"

void dasmerge( nchan, nsub, f_cen, f_inc, sp_in, 
              trunc, merge, clip,
              freq, sp_out, vel_out, numvals )
/*------------------------------------------------------------------------*/
/*
** Despike and merge (no DC adjustment) DAS spectra
*/
  int    nchan;          /* (in)   number of spectral channels (dim 'in') */
  int    nsub;           /* (in)   number of subbands                     */
  double f_cen[MAXBND];  /* (in)   observation (centre) frequency GHz     */
  double f_inc[MAXBND];  /* (in)   frequency step MHz                     */
  float  sp_in[MAXSPX];  /* (in)   observed spectrum                      */
  float  trunc;          /* (in)   number of channels/overlap to truncate */
                         /*        # = 0:     use defaults DTRUN1/2       */
                         /*        0 < # < 1: keep # fraction of overlap  */
                         /*        # >= 1:    drop # channels from ends   */
  int    merge;          /* (in)   0: no merge, else merge subbands       */
  float  clip;           /* (in)   amplitude cutoff (T<# or T>#: T=0)     */
                         /*        = 0 : use default DCLIP                */
  double freq[MAXSPX];   /* (out)  output frquencies     (GHz)            */
  float  sp_out[MAXSPX]; /* (out)  merged output spectrum                 */  
  float  vel_out[MAXSPX]; /* (out) output velocities     (km/s)           */
  int   *numvals;        /* (out)  Size of merged spectrum                */
/*------------------------------------------------------------------------*/
{
  int    i, j, k, l, m, n;
  int    nf,nval, nref;
  int    iweight[MAXSPX];
  float  overlap;
  double fcen, finc, fcen2, finc2, flo, fhi, fvel;
  double ddum;

  for (k = 0; k < MAXSPX; k++)
    {
      freq[k] = 0.0;
      sp_out[k]  = 0.0;
      vel_out[k] = 0.0;
      iweight[k] = 0.0;
    }


  nf = nchan;                         /* number of channels               */
  nval = nchan / nsub;                /* number of channels per subband   */   
  nref = nchan/2-1;                   /* vel. ref. channel                */

  if (clip < 0.00001)                 /* default CLIP                     */
    clip = DCLIP;

  /*
  ** Try determine the overlap from subband 1 and 2
  */
  if (nsub == 1)                       /* one subband only: trunc chan    */
    {
      overlap = 0.0;
      merge = 0;
      if (trunc < 1) trunc = DTRUN2;                     /* adopt default */
    }
  else
    {
      fcen  = f_cen[0];
      fcen2 = f_cen[1];
      finc  = (double) f_inc[0];
      ddum  =  ( (1000.0 * fcen2 - finc * (nval/2-1)) -
                 (1000.0 * fcen  + finc * (nval/2))  ) / finc;
      if (ddum < 0 ) ddum = -1.0 * ddum;
      overlap  = 2.0 * ( (int) ddum / 2.0 );
      if (trunc <= 0.001) trunc = DTRUN1;                /* adopt default */
    }

  if (trunc < 1) 
    trunc = 0.5 * (1.0-trunc) * overlap;                             
  else	
    if (trunc > (0.45*nval)) trunc = 0.45*nval;       /* don't trunc all! */

  trunc = (int) trunc;

  flo =  1.0e4;                      /* determine total span in frequency */
  fhi = -1.0e4;

  for (l = 0; l < nsub; l++)         /* cycle over subbands               */
  {
    ddum = 1.e-03*(1000.0 * f_cen[l] - (nval/2.0-1.0) * (double) f_inc[l]);
    if (ddum < flo) flo = ddum;
    if (ddum > fhi) fhi = ddum;    
    ddum = 1.e-03*(1000.0 * f_cen[l] + (nval/2.0) * (double) f_inc[l]);
    if (ddum < flo) flo = ddum;
    if (ddum > fhi) fhi = ddum;    
  }
   
  finc2 = (double) f_inc[(int) nsub/2];              /*  output increment */
  if (finc2 < 0) finc2 = -1.0 * finc2;               
  nf = nchan;                                       /* no output channels */
  if (merge != 0) 
    {
      nf = 1000.0 * (fhi-flo) /  finc2 + 1;
      if (nf > nchan) nf = nchan;
    }


  fvel = (f_cen[0]+f_cen[nsub-1])/2.0;             /* assume centred vel */


#ifdef DEBUG
  fprintf(stderr,"Number of bands: %d  Overlap: %f\n  Trunc: %f\n",
	  nsub, overlap, trunc);
  fprintf(stderr,"Nr output chans: %d  Merge: %d  Cut: %f\n",
	  nf, merge, clip);
  fprintf(stderr,"Flo: %f  Fhi: %f  Finc: %f\n",
          flo, fhi, finc2);
  fprintf(stderr,"Nref %4d  Fvel %f\n", nref, fvel);
#endif

  for (l = 0; l < nsub; l++)                       /* cycle over subbands */
    {
      j = l * nval;
      fcen = f_cen[l];
      finc = (double) f_inc[l]; 
      for (k = 0; k < nval; k++)                   /* cycle over channels */
	{
          ddum = 1.e-03 * (1000.0 * fcen + (k-nval/2.0+1.0) * finc);
	  if (k >= trunc && k < (nval - trunc))
	    if ( sp_in[j] > (-1.0 * clip) && sp_in[j] < (1.0 * clip) ) 
	      if (merge != 0)
		{
		  n = (1000.0 * (ddum-flo) / finc2 + 0.5);
		  if (n > (nf-1)) n = nf-1;
		  sp_out[n] += sp_in[j];
		  iweight[n] += 1;
		}
	      else
		{
		  sp_out[j] = sp_in[j];
		  iweight[j] = 1;
		}
          j++;
	}
    }


  for (n = 0; n < nf; n++) 
    if (iweight[n] != 0.0)
      sp_out[n] = sp_out[n]/iweight[n];

  if (merge != 0)                     /* make sure some variables (still) */
    nsub = 1;                         /* correct for multi and single     */
  else                                /* subband data                     */
    nf = nchan/nsub;

  for (l = 0; l < nsub; l++)
    {
      nval = nf;

#ifdef DEBUG
      fprintf(stderr," Write subband %1d data\n",l+1);
#endif

/*
** The central channel is exactly nchan / 2 because correlator data ran from
** chan '0' to chan 'nchan' (inclusive) but on writeout chan 0 got dropped.
** The frequency refers to the central frequency.
** 
*/

      fcen = f_cen[l];
      finc = (double) f_inc[l];

      for (k = 0; k < nval; k++)
	{
          m = l*nval+k;
	  if (merge != 0) 
	    freq[m] = flo + 1.0e-03 * m * finc2;
          else
	    freq[m] = 1.e-03 * (1000.0 * fcen + (k-nval/2.0+1.0) * finc);

          vel_out[m] = 2.997925e5 * (fvel - freq[m]) / fvel;

	}
    }

  *numvals = nf*nsub;
/* 
printf ("Returning out of dasmerge with nval = %d nsub = %d\n",nval, nsub);
printf ("Actvals in C prog = %d\n",*numvals); 
*/

}
