/*
  kperc.c            implements the kernel perceptron           

  2000 written by Ralf Herbrich                                 
  Technical University Berlin                                   

  2001 completed by Ralf Herbrich
  Microsoft Research Cambridge
   
  (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.
*/

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


/* macro for kernel calculations:
     a        pointer to vector of n components
     b        pointer to vector of n components
     k        counter variable
     n        number of dimensions
     K        inner product
     k_type   kernel type
     k_param  kernel parameter
*/
#define KERNEL(a,b,k,n,K,k_type,k_param) switch (k_type) { case 1: /* linear kernel */ for (K = .0, k = 0; k < n; k++) K += ((double) a [k] * (double) b [k]); break; case 2: /* polynomial kernel */ for (K = .0, k = 0; k < n; k++) K += ((double) a [k] * (double) b [k]); K = pow (K + 1.0, k_param); break; case 4: /* RBF kernel */ for (K = .0, k = 0; k < n; k++) K += (((double) a [k] - (double) b [k]) * ((double) a [k] - (double) b [k])); K = exp (-K/ (2.0 * k_param * k_param)); break; } 

/* implements one epoch of the kernel perceptron */
void kperc (int    *_n,       /* number of iterations */
	    int    *_idx,     /* current index to start at */
	    int    *_nochange,/* number of runs since last change */
	    int    *_l,       /* number of training points */
	    double *out,      /* cache of outputs */
	    double *alpha,    /* current state of the coefficients */
	    int    *y,        /* vector of classes (l x l) */
	    double *K,        /* kernel matrix (l x l) */
	    double *_eta,     /* learning rate */
	    double *_margin,  /* minimal required real valued output */
	    int    *_method,  /* learning method
			           0 - standard perceptron learning
				   1 - greedy standard perceptron learning
			      */
	    int    *mistakes, /* currrent number of mistakes */
	    int    *converged)/* indicates convergence */
{
  register int      i, iter;
  int               l = *_l, idx = *_idx, n = *_n;
  int               nochange = *_nochange, method = *_method;
  double            eta = *_eta, margin = *_margin, min, upd;


  /* in general, the routine has not converged */
  *converged = 0;
  
  /* M A I N   L O O P */
  for (iter = 0; iter < n; iter++) {
    
    if ((double) (y [idx]) * out [idx]  <= margin) {
                                              /* in case of misclassificaion */
      (*mistakes)++;                          /* count mistakes */

      /* standard learning rule */
      upd = eta * (double) (y [idx]);
      
      /* update alpha vector */
      alpha [idx] += upd; 
      /* update cached outputs */
      for (i = 0; i < l; i++)
	out [i] += upd * K [idx + i*l];
      nochange = 0;                           /* start counting again */
    } else {
      nochange++;                             /* increase counter */
    }

    if (method == 1) {
      for (idx = 0, min = (double) (y [0]) * out [0], i = 1; i < l; i++)
	if ((double) (y [i]) * out [i] < min) {
	  idx = i;
	  min = (double) (y [i]) * out [i];
	}
      printf ("idx = %d\n", idx);
    } else {      
      idx = (idx == (l-1))?0:idx+1;              /* set to next example */
    }

    if (nochange == l) {                       /* finish if there are no */
                                               /* more mistakes */
      *converged = 1;                          /* indicate convergence */
      break;                                   /* no more mistakes */
    }
  }

  /* return parameters */
  *_idx = idx;
  *_nochange = nochange;
  *_n = iter;
  
  return;
}

