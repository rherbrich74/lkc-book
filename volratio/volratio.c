/*
  volratio.c        numerical computation of the volume ratio

  2000 written by Ralf Herbrich
  Microsoft Research Cambridge
  
 (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.
 */

#include <stdio.h>
#include <math.h>

/*
  this function computes the real integral of \int_0^{acos(1-b)} (sin (x))^d dx
  by the rectangle method
*/
double nifty_integral (double b, int d, int no_rect) {
  double         res = 0.0;
  double         dx = (acos (1-b)) / (double) no_rect;
  double         x = dx / 2.0;
  register int   i;

  /* main loop */
  for (i = 0; i < no_rect; i++) {
    res += (pow (sin (x), (double) d)) * dx;
    x = x + dx;
  }

  /* return the integral value */
  return (res);
}

#define PREC      10000

/*
  interface to R 
 */
void volratio (double *b,             /* n x 1 vector of upper bounds */
	       int    *_d,            /* exponent on the sin */
	       int    *_n,            /* number of possible uppe bound values */
	       double *res)           /* n x 1 vector for the results */
{
  register int   i, j, k;        
  int            n = *_n, d = *_d;
  double         numerator, gamma;


  numerator = nifty_integral (2.0, d, PREC);
  for (i = 0; i < n; i++) {
    res [i] = log (numerator / nifty_integral (b [i], d, PREC));
  }
  
  return;
}


