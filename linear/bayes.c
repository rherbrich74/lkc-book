/*
   bayes.c            implements the Bayes Point Machines        
                                                               
   1999 written by Ralf Herbrich                                 
   Technical University Berlin

   2001 completed by Ralf Herbrich
   Microsoft Research Cambridge
   
   (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* computes the inner product of two vectors w.r.t the kernel K */
double get_inner (int l, double *alpha, double *beta, double *K)
{
  register int         i, j, k;
  double               inner;
  
  for (i = 0, k = 0, inner = .0; i < l; i++)
    for (j = 0; j < l; j++, k++)
      inner += (alpha [i] * beta [j] * K [k]);
  
  return (inner);
}

/* computes the SQUARED norm of weight vector alpha w.r.t. the kernel K */
double get_norm (int *_l, double *alpha, double *K, double *norm)
{
  return (*norm = get_inner (*_l, alpha, alpha, K));
}

/* normalize the weights to unity w.r.t. kernel K */
void normalize_weights (int    *_l,     /* number of training points */
			double *alpha,  /* vector to be normalized
					   (will be overwritten with the result) */
			double *K)      /* kernel matrix (l x l) */
{
  int                     l = *_l;
  register int            i;
  double                  dummy, norm = sqrt (get_norm (_l, alpha, K, &dummy));

  for (i = 0; i < l; i++) alpha [i] /= norm;
  return;
}

/* generates a random direction vector of unit length pointing
   TOWARDS the version space bounced at point m */
void random_version_direction (int    *_l,      /* number of training points */
			       int    *_m,      /* bounced at m */
			       double *K,       /* kernel matrix (l x l) */
			       int    *y,       /* vector of classes (l x 1) */
			       double *alpha)   /* resulting random vector (l x 1) */
{
  register int             i;
  int                      l = *_l, m = *_m, k;
  double                   inner;
  
  if (m < 0 || m > l)
    for (i = 0; i < l; i++) alpha [i] = 1.0 - 2.0 * drand48 ();
  else {
    /* repeat until a vector is found which points towards the version space */
    do {
      for (i = 0, k=m, inner = .0; i < l; i++, k+=l) {
	alpha [i] = 1.0 - 2.0 * drand48 ();
	inner += alpha [i] * (double) y [m] * K [k];
      }
    } while (inner <= .0);
  }
  
  normalize_weights (_l, alpha, K);
  return;
}

/* generates a random direction vector of unit length */
void random_direction (int    *_l,      /* number of training points */
		       double *K,       /* kernel matrix (l x l) */
		       int    *y,       /* vector of classes (l x 1) */
		       double *alpha)   /* resulting random vector (l x 1) */
{
  int         m = -1;
  
  random_version_direction (_l, &m, K, y, alpha);
  return;
}

/* computes the weighting coefficients in a "circular" geometry */
void compute_arc_avg (double frac, double inner, double *mu, double *lambda)
{
  *mu = -(frac*frac - frac*frac*inner-2.0)/(inner + 1.0);

  *mu = sqrt (*mu)*frac;
  *lambda = frac*frac*(1.0-inner) - 1.0;
  if (-(*mu)*inner + *lambda < .0)
    *lambda = -(*mu)*inner - *lambda;
  else
    *lambda = -(*mu)*inner + *lambda;    
  return;
}

#define TAU_MAX             10000         /* maximum flight time */

/* implements one step of the kernel billiard game */
void kbilliard (int    *_n,       /* number of iterations */
		int    *_l,       /* number of training points */
		double *alpha,    /* position of the ball (l x 1)
				     (will be overwritten with the result) */
		double *alpha2,   /* new position of the ball (l x 1) */
		double *alphad,   /* new distance vector (l x 1) */
		double *alphas,   /* sum of the two vectors (l x 1) */
		double *beta,     /* direction vector (l x 1) 
				     (will be overwritten with the result) */
		double *gamma,    /* estimate of the centre of mass (l x 1 
				     (will be overwritten with the result)) */
		double *K,        /* kernel matrix (l x l) */
		int    *y,        /* vector of classes (l x l) */
		double *L,        /* current length of trajectory */
		double *max_l,    /* maximum length of trajectory
				     (will be overwritten) */
		double *min_angle,/* minimum angle
				     (will be overwritten) */
		double *max_mu,   /* maximum mu (for abort criterion) */
		int    *last_min) /* index of the last hitted training point */
{
  double            add_l,         /* the length of current trajectory */
                    tau,           /* the actual flight time */
                    tau_min,       /* the minimal flight time */
                    d,             /* the distance of ball from hyperplane */
                    nu,            /* the projected direction vector */
                    nu_min = .0,   /* the minimal projected di... */
                    mu,            /* weighting coeff. of new ball position */
                    angle,         /* angular distance of centre & midpoint */
                    dummy,         /* dummy for get_norm */
                    lambda;        /* weighting coeff. of old centre of mass */
  register int      i, j, k, epoch;
  int               l = *_l, n = *_n, min;


  for (epoch = 0; epoch < n; epoch++) {
    
    /* compute next hyperplane for collision */
    do {
      tau_min = TAU_MAX;
      
      for (i = 0; i < l; i++) {
	for (j = 0, d = .0, nu = .0, k=i*l; j < l; j++, k++) {
	  d += (alpha [j] * K [k]);
	  nu += (beta [j] * K [k]);
	}
	d  *= ((double) y [i]);
	nu *= ((double) y [i]);
	
	/* consistency check! */
	if (d < -1e-4) {
	  fprintf (stderr, "d = %g by %d (last hit at %d)\n", d, i, *last_min);
	  d = .0;
	}
	tau = - d / nu;
	
	if ((tau > .0) && (tau < tau_min) && (i != *last_min)) {
	  tau_min = tau;
	  min = i;
	  nu_min = nu;
	}
      }
      
      /* generate a random direction vector towards the version space
	 in the case of escapes */
      if (tau_min == TAU_MAX) 
	random_version_direction (_l, last_min, K, y, beta);
      
    } while (tau_min == TAU_MAX);
    
    /* update the position of the ball */
    for (i = 0; i < l; i++) alpha2 [i] = alpha [i] + tau_min * beta [i];
    
    /* update the direction vector */
    beta [min] = beta [min] -  
      2.0 * (double) y [min] * nu_min / K [min * l + min];
    
    /* project the ball back onto the surface of the hypersphere */
    normalize_weights (_l, alpha2, K);
    
    /* make the direction vector of unit length */
    normalize_weights (_l, beta, K);
    
    /* update the center of mass of the trajectory */
    for (i = 0; i < l; i++) {
      alphad [i] = (alpha2 [i] - alpha [i]);
      alphas [i] = (alpha2 [i] + alpha [i]);
    }
    normalize_weights (_l, alphas, K);
    
    add_l = sqrt(get_norm (_l, alphad, K, &dummy));
    angle = get_inner (l, gamma, alphas, K);
    if (*max_l < add_l) *max_l = add_l;
    if (*min_angle > angle) *min_angle = angle;
    compute_arc_avg (*L/(*L + add_l), angle, &lambda, &mu);
    
    for (i = 0; i < l; i++) gamma [i] = lambda * gamma [i] + mu * alphas [i];
    
    /* update ball position and trajectory length */
    for (i = 0; i < l; i++) alpha [i] = alpha2 [i];
    *L = *L + add_l;
    *last_min = min;
    
    compute_arc_avg (*L / (*L + *max_l), *min_angle, &lambda, &mu);
    *max_mu = mu;
  }
  
  return;
}

