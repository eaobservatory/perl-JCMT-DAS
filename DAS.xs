/*        -*- C -*-

    JCMT data reduction subroutines for WORF2

	Done - dasmerge

						timj@jach.hawaii.edu

*/

#include "EXTERN.h"   /* std perl include */
#include "perl.h"     /* std perl include */
#include "XSUB.h"     /* XSUB include */

#include "arrays.h"
#include "dasmerge.h"  /* dasmerge header file */

MODULE = JCMT::DAS     PACKAGE = JCMT::DAS

void
dasmerge(nchan, nsub, f_cen, f_inc, sp_in, trunc, merge, clip, freq, sp_out, vel_out,numvals)
  int	nchan
  int	nsub
  double *	f_cen
  double * 	f_inc
  float *	sp_in
  float		trunc
  int 	merge
  float	clip
  double *	freq=NO_INIT
  float *	sp_out=NO_INIT
  float *	vel_out=NO_INIT
  int           numvals=NO_INIT
  CODE:
    sp_out = get_mortalspace(MAXSPX,'f');
    freq = get_mortalspace(MAXSPX,'d');
    vel_out = get_mortalspace(MAXSPX,'f');
    dasmerge(nchan, nsub, f_cen, f_inc, sp_in, trunc, merge, clip,
             freq, sp_out, vel_out, &numvals);
    unpack1D( (SV*)ST(8), (void *)freq, 'd', numvals);
    unpack1D( (SV*)ST(9), (void *)sp_out, 'f', numvals);
    unpack1D( (SV*)ST(10), (void *)vel_out, 'f', numvals);
  OUTPUT:
  numvals
