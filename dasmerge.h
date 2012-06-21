#define MAXBND    8
#define MAXSPX 2048

#define DTRUN1  0.5               /* Default fraction of overlap to keep */
#define DTRUN2   50               /* Default number of channels to trunc */
#define DCLIP  1000               /* Default clip level                  */

void dasmerge( int nchan, int nsub, double f_cen[], double f_inc[],
               float sp_in[], float trunc, int merge, float clip,
               double freq[], float sp_out[], float vel_out[], int * numvals );
