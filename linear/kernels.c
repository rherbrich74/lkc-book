/*
  kernels.c           Kernel routines for R

  1999 written by Ralf Herbrich
  Technical University Berlin

  2001 completed by Ralf Herbrich
  Microsoft Research Cambridge
   
  (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define K(a,X,Y,r,n,j,k,t,d,c,s)          a = .0; switch (type) { case 1: for (i = 0; i < n; i++) a += (X [i*r + j] * Y [k*n + i]); break; case 2: for (i = 0; i < n; i++) a += (X [i*r + j] * Y [k*n + i]); a = pow (a + 1.0, d); break; case 3: for (i = 0; i < n; i++) a += (X [i*r + j] * Y [k*n + i]); a = pow (a + c, d); break; case 4: for (i = 0; i < n; i++) a += ((X [i*r + j] - Y [k*n + i]) * (X [i*r + j] - Y [k*n + i])); a = exp (-(a)/(2.0 * s * s)); break; default: a = 0.0; break; }


/* this function is going to compute the inner product between a matrix
   and a matrix wrt. the choosen kernel */
void kinner (double *X,           /* a   r x n  matrix */
	     double *Y,           /* a   n x c  matrix*/
	     int    *_r,          /* number of rows of the X matrix */
	     int    *_n,          /* number of colums of the X matrix */ 
	     int    *_c,          /* number of colums of the Y matrix */
	     int    *_type,       /* type of the kernel involved */
	     double *_degree,     /* the degree of the polynomial */
	     double *_ofs,        /* the offset of the complete polynomial */
	     double *_sigma,      /* the bandwith of the RBF kernel */
	     int    *_norm,       /* normalise the kernel? */
	     double *ret)         /* return value */
{
  register int   i, j, k;
  int            r=*_r, c=*_c, n=*_n, type=*_type, norm=*_norm;
  double         degree=*_degree, sigma=*_sigma, ofs=*_ofs;
  double         *norm_X, *norm_Y, *tmp, *pt;
  
  /* if we need to normalise the kernel then we precompute the norm of
     each vector in the matrix X and Y */
  if (norm) {
    norm_X = (double *) malloc (r * sizeof (double));
    norm_Y = (double *) malloc (c * sizeof (double));
    tmp = (double *) malloc (n * sizeof (double));

    /* compute the norm of all the rows of X */
    for (j = 0; j < r; j++) {
      for (i = 0; i < n; i++) tmp [i] = X [i*r + j];
      K (norm_X [j], tmp, tmp, 1, n, 0, 0, type, degree, ofs, sigma);
      norm_X [j] = sqrt (norm_X [j]);
    }

    /* compute the norm of all the rows of Y */
    for (k = 0; k < c; k++) {
      pt = & (Y [k*n]);
      K (norm_Y [k], pt, pt, 1, n, 0, 0, type, degree, ofs, sigma);
      norm_Y [k] = sqrt (norm_Y [k]);
    }    
  }

  /* main loop for evalutation the kernel */
  for (j = 0; j < r; j++) {
    for (k = 0; k < c; k++) {
      K (ret [k*r + j], X, Y, r, n, j, k, type, degree, ofs, sigma);
      if ((norm) && (norm_X [j] != .0) && (norm_Y [k] != .0))
	ret [k*r + j] = ret [k*r + j] / (norm_X [j] * norm_Y [k]);
    }
  }

  /* free the memory */
  if (norm) {
    free (norm_X);
    free (norm_Y);
    free (tmp);
  }

  return;
}
