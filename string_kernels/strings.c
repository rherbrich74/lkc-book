/*
 * strings.c        program to demonstrate the ridge effect of
 *                  string kernels
 *
 * 2001 written by Ralf Herbrich
 * Microsoft Research Cambridge
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>

#define BUFFER_SIZE    102400

/******************************************************************/
/*  R E C U R S I V E     S U B S E Q U E N C E     K E R N E L   */
/******************************************************************/

double     ***cache;       /* this dynamic variable saves the auxillary
			      kernel values computed */

double Kprime (char *u, int p, char *v, int q, int n, double lambda) {
  register int         j;
  double               tmp;

  /* case 1: if a full substring length is processed, return*/
  if (n == 0) return (1.0);

  /* check, if the value was already computed in a previous computation */
  if (cache [n] [p] [q] != -1.0) return (cache [n] [p] [q]); 
  
  /* case 2: at least one substring is to short */
  if (p < n || q < n) return (0.0);
    
  /* case 3: recursion */
  for (j= 0, tmp = 0; j < q; j++) {
    if (v [j] == u [p - 1]) 
      tmp += Kprime (u, p - 1, v, j, n - 1, lambda) *   
        pow (lambda, (float) (q - j + 1));
  }

  cache [n] [p] [q] = lambda * Kprime (u, p - 1, v, q, n, lambda) + tmp;
  return (cache [n] [p] [q]);
}

double K (char *u, int p, char *v, int q, int n, double lambda) {
  register int  j;
  double        KP;

  /* the simple case: (at least) one string is to short */
  if (p < n || q < n) return (0.0);

  /* the recursion: use Kprime for the t'th substrings*/
  for (j = 0, KP = 0.0; j < q; j++) {
    if (v [j] == u [p - 1]) 
      KP += Kprime (u, p - 1, v, j, n - 1, lambda) * lambda * lambda;
  }
  
  return (K (u, p - 1, v, q, n, lambda) + KP);
}

/* recursively computes the subsequence kernel between s and t
   where subsequences of EXACTLY length n are considered */
double subsequence (char *u, char *v, int n, double lambda) {
  int           p = strlen (u), q = strlen (v), i, j, k;
  double        ret;

  /* allocate memory for auxiallary cache variable */
  cache  = (double ***) malloc (n * sizeof (double **));
  for (i = 1; i < n; i++) {
    cache  [i] = (double **) malloc (p * sizeof (double *));
    for (j = 0; j < p; j++) {
      cache  [i] [j] = (double *) malloc (q * sizeof (double));
      for (k = 0; k < q; k++) 
	cache  [i] [j] [k] = -1.0;
    }
  }
  
  /* invoke recursion */
  ret = K (u, p, v, q, n, lambda);
  
  /* free memory */
  for (i = 1; i < n; i++) {
    for (j = 0; j < p; j++) 
      free (cache  [i] [j]);
    free (cache  [i]);
  }
  free (cache);
  
  /* return the computed value */
  return (ret);
}

/******************************************************************/
/*                  S U B S T R I N G     K E R N E L             */
/******************************************************************/

/* efficiently compute the substring kernel between s and t
   where substrings up to length n are considered */
double full_substring (char *u, char *v, int n, double lambda) {
  int           p = strlen (u), q = strlen (v);
  register int  i, j, k;
  double        ret, tmp;

  /* computes the substring kernel */
  for (ret = 0.0, i = 0; i < p; i++) {
    for (j = 0; j < q; j++) 
      if (u [i] == v [j]) {
	for (k = 0, tmp = lambda * lambda;     /* starting condition */
	     (i + k < p) && (j + k < q) &&
	       (u [i + k] == v [j + k]) &&
	       (k < n);                        /* stop conditions */
	     k++, tmp *= (lambda * lambda))    /* update per iteration */
	  ret += tmp;
      }
  }
    
  /* return the computed value */
  return (ret);
}

/* efficiently compute the substring kernel between s and t
   where substrings of EXACTLY length n are considered */
double substring (char *u, char *v, int n, double lambda) {
  int           p = strlen (u), q = strlen (v);
  register int  i, j, k;
  double        ret, tmp;

  /* computes the substring kernel */
  for (ret = 0.0, i = 0; i < p; i++) {
    for (j = 0; j < q; j++) {
      for (k = 0, tmp = lambda * lambda;     /* starting condition */
	   (i + k < p) && (j + k < q) &&
	     (u [i + k] == v [j + k]) &&
	     (k < n);                        /* stop conditions */
	   k++, tmp *= (lambda * lambda));   /* update per iteration */
      
      if (k == n) ret += tmp;                /* update features in
						case of full match */
    }
  }
    
  /* return the computed value */
  return (ret);
}

/******************************************************************/
/*              B A G   O F   W O R D   K E R N E L               */
/******************************************************************/

/* efficiently compute the bag of word kernel between s and t
   where substrings of EXACTLY length n are considered */
double bow (char *u, char *v, int n, double lambda) {
  int           p = strlen (u), q = strlen (v);
  register int  i, j, k;
  double        ret, tmp;

  /* computes the bag-of-words kernel */  
  for (ret = 0.0, i = 0; i < p; i++) {
    for (j = 0; j < q; j++) {
      for (k = 0, tmp = lambda * lambda;     /* starting condition */
	   (i + k < p) && (j + k < q) &&     /* 1) there are still chars */
	     (u [i + k] == v [j + k]) &&     /* 2) they are equal */
	     (u [i + k] != ' ') &&           /* 3) they are unequal to ' ' */
	     (k < n);                        /* stop conditions */
	   k++, tmp *= (lambda * lambda));   /* update per iteration */

      if (k == n)                            /* substring detected */
	ret += tmp;
      for (;(j < q) && (v [j] != ' '); j++); /* skip the whole word in v */
    }

    for (;(i < p) && (u [i] != ' '); i++);   /* skip the whole word in u */
  }
    
  /* return the computed value */
  return (ret);
}

/*************************************************************/
/*               M A I N   R O U T I N E                     */
/*************************************************************/

#define UNDEFINED        -1
#define STRING           1
#define FULLSTRING       2
#define SEQUENCE         3
#define BOW              4

void string_kernel (char   **_base,    /* base directory of files */
		    int    *_no_docs,  /* number of documents to proceed */
		    char   **_type,    /* kernel type (as string) */
		    int    *_length,   /* substring length */
		    double *_lambda,   /* weight of each letter */
		    double *G) {       /* m x m Gram matrix (return) */
  int          depth=*_length;         /* depth of the string kernel */
  int          no_docs=*_no_docs;      /* no of documents */
  double       lambda=*_lambda;        /* decay factor */
  char         *base = *_base;         /* copy the base directory */
  char         *type = *_type;         /* copy the kernel type */
  int          k_type = UNDEFINED;     /* kernel type */
  
  int          i, j, k, found;      /* temporary files */
  int          f_size;              /* file size */
  FILE         *fp;                 /* temporary file pointer */
  char         name [256];          /* temporary file names */
  char         **docs;              /* array of document strings */
  char         tmp_s [BUFFER_SIZE]; /* temporary buffer for one document */


  /* check kernel type */
  if (strcmp (type, "substring") == 0) 
    k_type = STRING;
  else
    if (strcmp (type, "fullsubstring") == 0) 
      k_type = FULLSTRING;
    else
      if (strcmp (type, "subsequence") == 0)
	k_type = SEQUENCE;
      else 
	if (strcmp (type, "bow") == 0)
	k_type = BOW;
	else {
	  printf ("Second argument wrong. Stop program\n");
	  exit (-1);
	}
  
  /* read all documents into memory */
  printf ("[Reading documents]\n");
  docs = (char **) malloc (sizeof (char *) * no_docs);
  for (i = 0; i < no_docs; i++) {
    sprintf (name, "%s/%d.txt", base, i + 1);
    fp = fopen (name, "r");
    f_size = fread (tmp_s, 1, BUFFER_SIZE, fp);
    docs [i] = (char *) malloc (sizeof (char) * (f_size + 1));

    /*
      1) only copy the alphanumerical characters
      2) delete extra spaces
      3) convert everything into upper case
    */
    for (j = 0, k = 0, found = 0; j < f_size; j++) {
      if (!iscntrl ((int) tmp_s [j])) {
	if (isspace ((int) tmp_s [j])) {
	  if (!found) {
  	    docs [i] [k++] = toupper (tmp_s [j]); 
	    found = 1;
	  }
	} else {
 	  docs [i] [k++] = toupper (tmp_s [j]); 
	  found = 0;
	}
      }
      else
	if (!found)
	  docs [i] [k++] = ' ';
    }
    docs [i] [k-1] = '\0';    /* the last character is always a space (?) */
  }
    
  /* and compute the Gram matrix */
  printf ("[Computing the Gram matrix]\n");
  for (i = 0; i < no_docs; i++) {
    for (j = 0; j <= i; j++) {
      switch (k_type) {
      case STRING:
	G [i*no_docs + j] = substring (docs [i], docs [j], depth, lambda);
	break;
      case FULLSTRING:
	G [i*no_docs + j] = full_substring (docs [i], docs [j], depth, lambda);
	break;
      case SEQUENCE:
	G [i*no_docs + j] = subsequence (docs [i], docs [j], depth, lambda);
	break;
      case BOW:
	G [i*no_docs + j] = bow (docs [i], docs [j], depth, lambda);
	break;
      }

      /* copy the value since the Gram matrix is symmetric */
      G [j*no_docs + i] = G [i*no_docs + j];
	
      printf ("\rProgress: %5.2f%%", (float) i * 100.0/ (float) no_docs);
      fflush (stdout);
    }
  }
  printf ("\n[Gram matrix computed]\n");

  /* free all memory */
  for (i = 0; i < no_docs; i++) 
    free (docs [i]);
  free (docs);
  
  return;
}
