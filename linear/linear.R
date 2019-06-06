### linear.R    a library for learning linear classifiers
###
### 1999/2000 written by Ralf Herbrich
### Technical University Berlin
###
### 2001 completed by Ralf Herbrich
### Microsoft Research Cambridge
###   
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
##     L O A D   E X T E R N A L   L I B R A R I E S
############################################################

source ("linear/kernels.R");

if (!is.loaded ("pr_loqo")) {
  cat ("Loading dynamic library for LOQO optimizer\n")
  if (R.version$os == "Win32") {
    dyn.load ("linear/pr_loqo.dll");
  } else {
    dyn.load ("linear/pr_loqo.so");
  }
}

if (!is.loaded ("kbilliard")) {
  cat ("Loading dynamic library for Bayes Point Machines\n")
  if (R.version$os == "Win32") {
    dyn.load ("linear/bayes.dll");
  } else {
    dyn.load ("linear/bayes.so");
  }
}

if (!is.loaded ("kperc")) {
  cat ("Loading dynamic library for Kernel Perceptrons\n")
  if (R.version$os == "Win32") {
    dyn.load ("linear/kperc.dll");
  } else {
    dyn.load ("linear/kperc.so");
  }
}

############################################################
##     P E R C E P T R O N   R O U T I N E S
############################################################

############################################################
## this routine defines a new perceptron
############################################################

perc.new <- function (X, y, k) {
  if (!is.matrix (X))
    stop ("X (1st arg) has to be a matrix");
  if (!is.vector (y))
    stop ("y (2nd arg) has to be a vector");
  if (dim (X) [1] != length (y))
    stop ("no of rows of X has to be length of y");
  if (class (k) != "kernel")
    stop ("k (3rd arg) has to be of class kernel");
  
  ret <- list (data=X,
	       class=y,
	       kernel=k,
	       K=kernel.inner (X, t(X), k),
	       alpha=rep (0, length(y)),
	       theta=0,
	       margin=0,
	       eps=1e-8,
               classifier=prod(abs (y)==1));
  class (ret) <- "perc";
  return (ret);
}

############################################################
## classifies a given set of points
############################################################

perc.classify <- function (X, p) {

  if (!is.matrix (X))
    stop ("X (1st arg) has to be a matrix");
  if (class (p) != "perc")
    stop ("p (2nd arg) has to be a perceptron");

  H <- kernel.inner (p$data, t(X), p$kernel);
  return (p$alpha %*% H + p$theta);
}

############################################################
## computes the classification error given a test set
############################################################

perc.class.error <- function (X, y, p) {
  if (!is.matrix (X)) 
    stop ("X (1st arg) has to be a matrix");
  if (!is.vector (y))
    stop ("y (2nd arg) has to be a vector");
  if (dim (X)[1] != length (y))
    stop ("number of rows of X and elements in y have to be equal");
  if (class (p) != "perc")
    stop ("p (3rd arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("classification error is meaningless for regression");    

  return ((sum (sign (perc.classify (X, p)) != y)) / length (y));
}

############################################################
## computes the sparsity in terms of the expansion coefficients
############################################################

perc.sparsity <- function (p) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");

  return (sum ((abs (p$alpha) > p$eps)));
}

############################################################
## computes the (classical) margin
############################################################

perc.margin <- function (p) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");

  ## normalise the weights to unit length
  norm <- sqrt (t(p$alpha) %*% p$K %*% p$alpha);
  out <- min (p$class * (p$K%*%p$alpha+p$theta));
  return (out / norm);
}

############################################################
## computes the (normalised) margin
############################################################

perc.new.margin <- function (p) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");

  ## normalise the weights AND data points in feature space
  ## to unit length
  norm <- sqrt (t(p$alpha)%*%p$K%*%p$alpha);
  out <- min (p$class * (p$K%*%p$alpha+p$theta) / sqrt (diag(p$K)));
  return (out / norm);
}

############################################################
## penalises the Gram matrix diagonal by lambda>0
############################################################

perc.incdiag <- function (p, lambda=0) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");

  diag (p$K) <- diag (p$K) + lambda;
  return (p);
}

############################################################
## prints all information about the perceptron onto stdout
############################################################

print.perc <- function (p, ...) {
  cat ("Perceptron\n=======\n\nAlpha's\n-------\n\n", ...);
  for (i in 1:length (p$alpha)) 
    cat ("[", i, "]\t", p$alpha[i],  "\n", ...);
  cat ("\n\nTheta      = ", p$theta,     "\n", ...);
  cat ("Margin     = ",     p$margin,    "\n", ...);
  cat ("EPS        = ",     p$eps,       "\n", ...);
  cat ("CLASSIFIER = ",     p$classifier,"\n", ...);
}

######################################################################
##  K E R N E L   P E R C E P T R O N   L E A R N I N G
######################################################################

perc.learn <- function (p, eta=1, n=1000000000,
                        method="standard", margin=0,
                        verbosity=0) {
  
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.learn can only be applied to classifiers");
  if ((method != "standard") && (method != "greedy"))
    stop ("perc.learn must use 'standard' or 'greedy' for method");
  if (margin < 0)
    stop ("perc.learn must be provided with positive margins");

  l <- length (p$alpha);                      # number of training points
  idx <- 0;                                   # next training point to
                                              # be considered
  nochange <- 0;                              # number of steps since
                                              # last change
  mistakes <- 0;                              # number of mistakes
  converged <- 0;                             # convergence indicator
  method.id <- switch (method,                # convert method "string" into
                                              # method "id"
                       "standard"=0,
                       "greedy"=1
                       );
  
  if (verbosity == 0) {
    res <- .C ("kperc",
               n=as.integer (n),
               idx=as.integer (idx),
               nochange=as.integer (nochange),
               as.integer (l),
               out=as.double (vector ("numeric", l)),
               alpha=as.double  (p$alpha),
               as.integer (p$class),
               as.double  (p$K),
               as.double  (eta),
               as.double  (margin),
               as.integer (method.id),
               mistakes=as.integer (mistakes),
               converged=as.integer (converged));    
    p$alpha <- res$alpha;
    p$margin <- min (p$class * (p$K %*% p$alpha));
    if (verbosity > 0) {
      if (res$converged == FALSE || p$margin < 0) {
        cat ("perceptron learning is not converged");
        browser ();
      }
    }
  }

  # watch the learning 
  if (verbosity == 1) {

    cat ("Invoke perceptron learning\n");
    res <- .C ("kperc",
               n=as.integer (n),
               idx=as.integer (idx),
               nochange=as.integer (nochange),
               as.integer (l),
               out=as.double (vector ("numeric", l)),
               alpha=as.double  (p$alpha),
               as.integer (p$class),
               as.double  (p$K),
               as.double  (eta),
               as.double  (margin),
               as.integer (method.id),
               mistakes=as.integer (mistakes),
               converged=as.integer (converged));    
    p$alpha <- res$alpha;
    p$margin <- min (p$class * (p$K %*% p$alpha));
    
    cat ("Finished perceptron learning\n\n");
    cat ("Steps     ", res$n, "\n");
    cat ("Mistakes: ", res$mistakes, "\n");
    cat ("Sparsity: ", perc.sparsity (p), "\n");
    cat ("Margin:   ", perc.new.margin (p), "\n");
  }
  
  return (p);
}

############################################################
##     S U P P O R T   V E C T O R   L E A R N I N G
############################################################

############################################################
## performs linear soft margin SVM's (if C=1e20 it virtually
## performs hard margin classification)
############################################################

perc.svm <- function (p,                # the original perceptron
                      r=0,              # the vector of the r.h.s.
                      eta=1,            # the vector of factors on l.h.s.
                      C=1000,           # trade-off constant
                      bias=FALSE,       # constant term included?
                      norm=FALSE,       # normalise the data in feature space?
                      verbosity=0) {  
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.svm can only be applied to classifiers");
  
  n <- length (p$alpha);                # number of variables to opt.
  
  if (norm == TRUE) {
    S <- diag (p$class * eta * 1/sqrt (diag (p$K)));
    H <- S %*% p$K %*% S;               # Hessian matrix
                                        # with normalisation!!!
  } else {
    H <- diag (p$class * eta) %*% p$K %*% diag (p$class * eta);
                                        # Hessian matrix
  }
  if (r == 0) {
    cc <- rep (-1,n);                   # linear ojective function
  } else {
    cc <- -r;                           # user defined importance vector
  }
  l <- rep (0,n);                       # lower bounds on variables
  u <- rep (C,n);                       # upper bounds on variables

  if (!bias) {
    m <- 0;                             # no constraints
    A <- 0;                             # dito
    b <- 0;                             # dito
  }
  else {
    m <- 1;                             # one constraint
    A <- p$class;                       # \sum_i y_i \alpha_i = 0
    b <- 0;
  }

  primal <- rep(0, 3*n);                # variables of the primal 
  dual   <- rep(0, m+2*n);              # variables of the dual

  p$alpha <- .C ("pr_loqo_R",
		 as.integer (n),
		 as.integer (m),
		 as.double  (cc),
		 as.double  (H),
		 as.double  (A),
		 as.double  (b),
		 as.double  (l),
		 as.double  (u),
		 primal=as.double (primal),
		 as.double  (dual),
		 as.integer (verbosity),
		 as.double  (7),
		 as.integer (10000),
		 as.double  (0.01),
		 as.double  (10),
		 as.integer (1))$primal [1:n] * p$class * eta;

  valid <- (C - abs (p$alpha)) > 1e-10; # determine correctly classified points
  if (norm ==TRUE) {
    p$alpha <- p$alpha / sqrt (diag (p$K));
  }
  
  if (bias) {
    p$theta <- 0;                            
    pos <- (perc.classify (p$data [p$alpha > 1e-10 & valid,], p));
    neg <- (perc.classify (p$data [p$alpha < -1e-10 & valid,], p));
    n <- min (length (pos), length (neg));
    pos <- mean (pos [1:n]); 
    neg <- mean (neg [1:n]);  
    p$theta <- -(pos + neg) / 2;        # theta = 1/2(mean(B) - mean(A))
  }
  else 
    p$theta <- .0;                      # else theta = .0

  alpha <- p$alpha [valid];
  p$margin <- min (r);                  # margin was fixed to 1
  p$eps <- min (abs (alpha [perc.classify (p$data [valid,], p) *
                             p$class [valid] - p$margin < 1e-4])); 
  

  return (p);
}

############################################################
## performs \nu-SVM learning
############################################################

perc.nu.svm <- function (p, nu=0.2, bias=FALSE, verbosity=0) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.nu.svm can only be applied to classifiers");

  n <- length (p$alpha);                 # number of variables to opt.
  H <- diag (p$class) %*% p$K %*% diag (p$class);    # Hessian matrix
  C <- rep (0,n);                        # linear objective function
  lower <- rep (0, n);                   # lower bounds
  upper <- rep (1/n, n);                 # upper bounds
  
  if (!bias) {                              
    m <- 1;                              # one constraint
    A <- rep (1, n);                     # \sum_i \alpha_i = \nu
    b <- nu;                             # see above
  }
  else {
    m <- 2;                              # two constraints
    A <- rbind (rep (1, n), p$class);    # \sum_i \alpha_i = \nu
    b <- c (nu, 0);                      # \sum_i y_i \alpha_i = 0
  }
  
  primal <- rep (0, 3*n);                # variables of the primal
  dual   <- rep (0, m+2*n);              # variables of the dual 

  p$alpha <- .C ("pr_loqo_R",
		 as.integer (n),
		 as.integer (m),
		 as.double  (C),
		 as.double  (H),
		 as.double  (A),
		 as.double  (b),
		 as.double  (lower),
		 as.double  (upper),
		 primal=as.double (primal),
		 as.double  (dual),
		 as.integer (verbosity),
		 as.double  (10),
		 as.integer (10000),
		 as.double  (0.01),
		 as.double  (10),
		 as.integer (1))$primal [1:n] * p$class;
  
  p$alpha <- p$alpha / max (abs (p$alpha)); # rescale alpha's for eps-determ.
  valid <- (1 - abs (p$alpha)) > 1e-10;     # determine non-errors
  p$theta <- 0;               
  pos <- (perc.classify (p$data [p$alpha > 1e-3 & valid,], p));
  neg <- (perc.classify (p$data [p$alpha < -1e-3 & valid,], p));
  n <- min (length (pos), length (neg));
  pos <- mean (pos [1:n]);                  # determine mean of positive/
  neg <- mean (neg [1:n]);                  # negative vectors
  
  if (bias) {
    p$theta <- -(pos + neg) / 2;            # theta = 1/2 (B-A)
  }
  else 
    p$theta <- .0;                          # no bias
  p$margin <- (pos - neg) / 2;              # margin = 1/2 (A-B) (adaptive!)
  
  alpha <- p$alpha [valid]; 
  p$eps <- min (abs (alpha [perc.classify (p$data [valid,], p) *
                            p$class[valid] - p$margin < 1e-4]));
                                            # determine accuracy via valid
  return (p);
}

######################################################################
##  G A U S S I A N   P R O C E S S E S
######################################################################

## regression estimation with Gaussian processes
perc.GP.regress <- function (p, var=0) {
  
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == TRUE)
    stop ("perc.GP.regress can only be applied to regression");
  if (var < 0)
    stop ("Variance (2nd arg) has to be positive");

  p$alpha <- solve (p$K + diag (rep (var, length=length (p$alpha))), p$class);
  p$margin <- 0;
  
  return (p);
}

## INTERNAL FUNCTION: never call directly!
## returns the sigmoidal transformation of t
sigmoid <- function (t, C=1) {
  return (1 / (1 + exp (-C * t)));
}

## INTERNAL FUNCTION: never call directly!
## returns the objective function value
GP.J <- function (alpha, y, out, G) {
  return (1/2 * t (y + 1) %*% out - sum (log (1 + exp (out))) -
          1/2 * t (alpha) %*% G %*% alpha)
}

## classification with Gaussian processes using the Laplace
## approximation
perc.GP.class <- function (p,            # old linear classifier
                           C=1,          # steepness of the sigmoid
                           var=0,        # variance on the
                                         # latent variable
                           verbosity=0,  # output convergence msgs?
                           tol=1e-8,     # stopping criterion
                           backoffs=8) { # 2^-backoffs is a
                                         # lower bound on eta
  
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.GP.class can only be applied to regression");
  if (var < 0)
    stop ("Variance (3rd arg) has to be positive");
  if (C <= 0)
    stop ("Steepness parameter must be strictly positive");
  if (tol <= 0)
    stop ("Tolerance must be strictly positive");
  if (backoffs <= 0)
    stop ("At least one backoff step must be allowed");

  ## extract training set size
  m <- length (p$alpha);
  
  ## create the vector and diagonal matrix of ones
  one <- cbind (rep (1, length=m));
  I <- diag (rep (1, length=m));

  ## copy the classes
  y <- cbind (p$class);
  
  ## compute Gram matrix (plus variance)
  G <- C * (p$K + diag (rep (var, length=m)));

  ## states of the latent hidden variables
  p$alpha <- matrix (0, nrow=m, ncol=1);
  out <- G %*% p$alpha;
  
  ## compute the function value
  J <- GP.J (p$alpha, y, out, G);
  if (verbosity) {
    cat ("J =", formatC (J), "\n");
  }

  ## NEWTON - RAPHSON optimisation
  while (TRUE) {
    ## compute the sigmoidial transformation of the
    ## current latent variables
    v <- cbind (sigmoid (out));
    
    ## compute the Gradient
    grad <- 1/2 * (y + one) - v - p$alpha;

    ## compute the Hessian
    P <- diag (as.vector (v * (1 - v)));
    H <- -(P %*% G + I);

    ## compute the offset vector
    delta <- solve (H, grad);

    ## adjust the stepsize of the offset vector
    eta <- 1;
    for (i in 1:backoffs) {
      ## do the update by eta
      alpha.tmp <- p$alpha - eta * delta;
      out.tmp <- G %*% alpha.tmp;
      J.tmp <- GP.J (alpha.tmp, y, out.tmp, G);

      ## check if the function has been increased
      if (J.tmp > J) {
        if (verbosity) {
          cat ("J =", J.tmp, "\t",
               sqrt (t (grad) %*% grad) / m, "\n");
        }
        
        ## update parameters
        J <- J.tmp;
        out <- out.tmp;
        p$alpha <- alpha.tmp;
        break;
      } else {
        eta <- eta / 2;
      }
    }
    
    if (i == backoffs) {
      cat ("Stopping due to back-off limit, gradient = ",
           sqrt (t (grad) %*% grad)/m, "\n");
      break;
    }

    ## stopping criterion
    if (sqrt (t (grad) %*% grad) / m < tol) {
      break;
    }
  }

  ## type conversion only
  p$alpha <- as.vector (p$alpha);
  p$margin <- 0;
  
  return (p);
}

######################################################################
##  R E L E V A N C E   V E C T O R     M A C H I N E S
######################################################################

## regression estimation with RVMs
perc.RVM.regress <- function (p,              # current perceptron
                              tol=sqrt(.Machine$double.eps), # tolerance used to
                                              # determine pruning
                              theta=1,        # start value of all theta's
                              var=0.1,        # variance of outputs
                              var.fix=FALSE,  # fixed variance?
                              iterations=100, # no. of iterations
                              verbosity=0) {  # log info.?

  ## parameter checks
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == TRUE)
    stop ("perc.RVM.regress can only be applied to regression");

  ## initialise variables
  m <- length (p$alpha);
  theta <- rep (theta, length=m);
  alpha <- rep (0, length=m);

  ## copy Gram matrix
  G <- p$K;

  ## precompute the inner products of the data matrix
  ## with the targets
  Xt <- t (G) %*% p$class;

  ## MAIN loop 
  for (i in 1:iterations) {
    ## determine which theta's are below machine precision
    non.zero <- (theta > tol);

    ## prune the theta's by equating them to zero and
    ## set the corresponding alpha's to zero as well
    theta [!non.zero] <- 0;
    alpha [!non.zero] <- 0;
    
    ## count the number of features used
    n <- sum (non.zero);

    ## reduce data matrix and coefficient vector
    X <- G [1:m, non.zero];
    theta.non.zero <- theta [non.zero];

    ## compute the Cholesky decomposition of the Hessian
    ## (\Sigma^{-1} in the book)
    ## (much more stable - according to Mike Tipping)
    Sigma.inv <- 1/var * t (X) %*% X + diag (1/theta.non.zero);
    R <- chol (Sigma.inv);
    
    ## invert the Cholesky decomposition
    R.inv <- solve (R);

    ## compute the new alpha coefficients
    alpha [non.zero] <- 1/var * (R.inv %*% (t (R.inv) %*% Xt [non.zero]))

    ## compute the diagonal of Sigma quickly
    diag.Sigma <- apply (R.inv^2, 1, sum);

    ## compute the error term
    error <- sum ((p$class - X %*% alpha [non.zero])^2);

    ## log some information
    if (verbosity > 0) {
      ## compute the log determinant of Sigma.inv using the fact
      ## that the Cholesky decomposition R is upper triangular matrix
      ## and thus it is just the sum of the log of the diagonal
      log.det.Sigma.inv <- 2 * sum (log (diag (R)));
      
      ## compute the evidence to monitor convergence
      evidence <- -1/2 * (log.det.Sigma.inv +
                          sum (log (theta.non.zero)) +
                          m * log (var) + 1/var * error +
                          (alpha [non.zero]^2) %*% (1/theta.non.zero));

      cat ("Evidence =", formatC (evidence), "\tn=", n, "\tvar=", var, "\n");

      if (verbosity > 1) {
        p$alpha <- alpha;
        plot (p, mag=0.0, grid=200);
      }
    }

    ## compute the auxillary variables \zeta
    zeta <- 1 - diag.Sigma / theta.non.zero;

    ## update the theta's
    theta [non.zero] <- alpha [non.zero]^2 / zeta;

    ## update the variance
    if (var.fix == FALSE) {
      var <- error / (m - sum (zeta));
    }
  }

  ## copy the soultion found
  p$alpha <- alpha;
  p$margin <- 0;
  
  return (p);
}

## INTERNAL FUNCTION: never call directly!
## returns the objective function value
RVM.J <- function (alpha, y, out, Theta.inv) {
  return (1/2 * t (y + 1) %*% out - sum (log (1 + exp (out))) -
          1/2 * t (alpha) %*% Theta.inv %*% alpha);
}

## classification with RVMs
perc.RVM.class <- function (p,              # current perceptron
                            tol=sqrt(.Machine$double.eps), # tolerance used to
                                            # determine pruning
                            theta=1,        # start value of all
                                            # theta's
                            iterations=100, # no. of iterations
                            C=1,            # steepness of the sigmoid
                            backoffs=8,     # 2^-backoffs is a lower bound
                                            # on eta
                            verbosity=0) {  # log info.?

  ## parameter checks
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.RVM.class can only be applied to classification");
  if (C <= 0)
    stop ("Steepness parameter must be strictly positive");
  if (tol <= 0)
    stop ("Tolerance must be strictly positive");
  if (backoffs <= 0)
    stop ("At least one backoff step must be allowed");

  ## initialise variables
  m <- length (p$alpha);
  theta <- rep (theta, length=m);
  alpha <- rep (0, length=m);
  one <- cbind (rep (1, length=m));
  y <- cbind (p$class);
  
  ## copy Gram matrix as the data description
  ## also incorporate the steepness directly here
  G <- C * p$K;

  ## MAIN loop 
  for (i in 1:iterations) {
    ############################################################
    ## 1. PRUNE THE (ALMOST) ZERO COEFFICIENTS
    ############################################################

    ## determine which theta's are below machine precision
    non.zero <- (theta > tol);

    ## prune the theta's by equating them to zero and
    ## set the corresponding alpha's to zero as well
    theta [!non.zero] <- 0;
    alpha [!non.zero] <- 0;
    
    ## count the number of features used
    n <- sum (non.zero);

    ## reduce data matrix and coefficient vector
    X <- G [1:m, non.zero];
    theta.non.zero <- theta [non.zero];

    ############################################################
    ## 2. COMPUTE THE LAPLACE APPROXIMATION (Newton-Raphson)
    ############################################################

    ## initialise the Newton-Raphson iteration
    out <- X %*% alpha [non.zero];
    Theta.inv <- diag (1/theta.non.zero);
    J <- RVM.J (alpha [non.zero], y, out, Theta.inv);
    
    ## log some information
    if (verbosity > 1) {
      cat ("[Invoking Laplacian approximation]\nJ =", formatC (J),"\n");
    }
    while (TRUE) {
      ## compute the sigmoidial transformation of the
      ## current latent variables
      v <- cbind (sigmoid (out));
    
      ## compute the Gradient
      grad <- t (X) %*% ((y + one) / 2 - v) - Theta.inv %*% alpha [non.zero];
      
      ## compute the Hessian
      P <- diag (as.vector (v * (1 - v)));
      H <- -(t (X) %*% P %*% X + Theta.inv);
      
      ## compute the offset vector
      delta <- solve (H, grad);

      ## adjust the stepsize of the offset vector
      eta <- 1;
      for (i in 1:backoffs) {
        ## do the update by eta (but copy alpha first)
        alpha.tmp <- alpha;
        alpha.tmp [non.zero] <- alpha [non.zero] - eta * delta;
        out.tmp <- X %*% alpha.tmp [non.zero];
        J.tmp <- RVM.J (alpha.tmp [non.zero], y, out.tmp, Theta.inv);

        ## check if the function has been increased
        if (J.tmp > J) {
          if (verbosity > 1) {
            cat ("J =", J.tmp, "\t",
                 sqrt (t (grad) %*% grad) / n, "\n");
          }
        
          ## update parameters
          J <- J.tmp;
          out <- out.tmp;
          alpha <- alpha.tmp;
          break;
        } else {
          eta <- eta / 2;
        }
      }
    
      if (i == backoffs) {
        if (verbosity > 0) {
          cat ("Stopping due to back-off limit, gradient = ",
               sqrt (t (grad) %*% grad)/n, "\n");
        }
        break;
      }

      ## stopping criterion
      if (sqrt (t (grad) %*% grad) / m < tol) {
        break;
      }
    }

    ############################################################
    ## 3. MAKE ONE UPDATE STEP 
    ############################################################

    Sigma.inv <- -H;
    R <- chol (Sigma.inv);
    
    ## invert the Cholesky decomposition
    R.inv <- solve (R);

    ## compute the diagonal of Sigma quickly
    diag.Sigma <- apply (R.inv^2, 1, sum);

    ## compute the auxillary variables \zeta
    zeta <- 1 - diag.Sigma / theta.non.zero;

    ## update the theta's
    theta [non.zero] <- alpha [non.zero]^2 / zeta;

    ## log some information
    if (verbosity > 0) {
      cat ("L =", formatC (J), "\tn=", n, "\n");
      if (verbosity > 1) {
        p$alpha <- as.vector (alpha);
        plot (p);
      }
    }
  }
  
  ## copy the soultion found
  p$alpha <- as.vector (alpha);
  p$margin <- 0;
  
  return (p);
}

############################################################
##     B A Y E S   P O I N T   M A C H I N E S
############################################################

perc.bayes  <- function (p,                # input perceptron
					   # (HAS TO BE IN VERSION SPACE!)
			 max.epoch = Inf,  # number of maximum epochs
			 tol = 1e-3,       # tolerance value of maximum update
			 disp.epoch = 100, # number of epochs for output (verb > 0)
                         test.X,           # a matrix of test observations
                         test.y,           # a vector of test classes
                         plot = FALSE,     # should we generate a plot
                         plot.svm = plot,  # should we generate a SVM circle plot?
			 verbosity = 0,    # verbosity level
                         ...) 
{
  if (p$theta != 0) 
    stop ("only perceptrons WITHOUT theta allowed");
  if (plot && (p$kernel$type != 1 || dim (p$data) [2] != 3)) {
    warning ("Turn off plot option (wrong kernel or wrong dimension of data)");
    plot <- FALSE;
  }
  if (verbosity == 2 && missing (test.X)) {
    warning ("Switch back to verbosity level 1 (no test data available)");
    verbosity <- 1;
  }
  if (p$classifier == FALSE)
    stop ("perc.bayes can only be applied to classifiers");
      
  l         <- length (p$alpha);        # number of training points
  ball      <- p$alpha;                 # the position of the ball
  balln     <- rep (0, length (ball));  # a temporary memory location for the 
					# new position of the ball
  balld     <- rep (0, length (ball));  # the difference vector of the ball
  balls     <- rep (0, length (ball));  # the sum vector of the ball
  dir       <- rep (0, length (ball));  # the direction vector
  p$alpha   <- rep (0, length (ball));  # the centre of mass (p$alpha)
  L         <- 0;                       # current length of trajectory
  max.l     <- 0;                       # maximum length of trajectory
  min.angle <- 1.0;                     # minimum angle (between 0 and 1)
  max.mu    <- 0.0;                     # maximum weighting factor
  last.min  <- -1;                      # index of last hitted training example
  stop.epoch<- disp.epoch;              # copy disp.epoch
  
  ## log some information
  switch (verbosity,
          cat ("\tEpoch\tTOL\t\tLast.min\n==========================================\n\n"),
          cat ("\tEpoch\ttest\tTOL\t\tLast.min\n=================================================\n\n"));

  ## initalize with random direction
  dir <- .C ("random_direction",
	     as.integer (l),
	     as.double  (p$K),
	     as.integer (p$class),
	     out = as.double  (dir))$out;

  ## normalize the weights to unit length
  ball <- .C ("normalize_weights",
	      as.integer (l),
	      out = as.double  (ball),
	      as.double  (p$K))$out;

  ##  in the case of plot, generate a new output
  VT <- list();
  track.ball <- 0;
  pt <- 0;
  svm.pt <- 0;
  if (plot) {
    svms <- abs (ball) > 1e-4;
    stop.epoch <- 1;
    plot.new ();
    track.ball <- t (t (p$data) %*% ball);
    oldpar <- par();
    VT <- plot.3d (ball (ntheta=30, nphi=30), pch=18, col="black",
                   axes=F, xlab="", ylab="", ...);
    svm.pt <- project.3d (track.ball, VT=VT);
    pt <- svm.pt;
  }
  
  epoch <- 0;
  while (TRUE) {
    ## print out some log information
    if (epoch %% disp.epoch == 0) {
      switch (verbosity,
              cat ("\t", epoch,
                   "\t", formatC (max.mu, format="f"),
                   "\t", last.min, "\n"),
              cat ("\t", epoch,
                   "\t", formatC (perc.class.error (test.X, test.y, p),
                                  format="f"),
                   "\t", formatC (max.mu, format="f"),
                   "\t", last.min, "\n"));
    }

    ## do a set of steps
    res <- .C ("kbilliard",
               as.integer (stop.epoch),
	       as.integer (l),
	       alpha   = as.double  (ball),
	       as.double  (balln),
	       as.double  (balld),
	       as.double  (balls),
	       dir    = as.double (dir),
	       centre = as.double (p$alpha),
	       as.double  (p$K),
	       as.integer (p$class),
	       L = as.double (L),
	       maxl = as.double (max.l),
	       mina = as.double (min.angle),
	       maxmu = as.double (max.mu),
	       lastmin = as.integer (last.min));
   
    ## copy the results into current vectors
    ball      <- res$alpha;
    dir       <- res$dir;
    p$alpha   <- res$centre;
    L         <- res$L;
    max.l     <- res$maxl;
    min.angle <- res$mina;
    max.mu    <- res$maxmu;
    last.min  <- res$lastmin;

    ## in the case of plot, print the output
    if (plot) {
      oldpt <- pt;
      pt <- t (t (p$data) %*% ball);
      track.ball <- rbind (track.ball, pt);
      pt <- project.3d (pt, VT=VT);
      points (pt [1], pt [2], pch=".", col="red");
      lines (c (oldpt [1], pt [1]),
             c (oldpt [2], pt [2]),
             col=par()$bg, lwd=2);
    }
    
    epoch <- epoch + stop.epoch;
    if (epoch >= max.epoch || max.mu < tol) {
      break;
    }
  }
  p$eps <- 0;
  p$margin <- 0;
  p$theta <- 0;

  ## if plot, the log the ball positions in "ball.dat"
  if (plot) {
    points (project.3d (track.ball [2:epoch,], VT=VT), col="red", pch=19);
    sol <- t (t (p$data) %*% p$alpha);
    points (project.3d (sol, VT=VT), col="black", pch=5, lwd=2, cex=2);
    points (svm.pt [1], svm.pt [2], col="black", pch=4, cex=2, lwd=2);
    save (track.ball, file="ball.dat", ascii=TRUE);

    ## if we also should plot the function ball of SVM's?
    w <- track.ball [1,];
    if (plot.svm) {
      gamma <- max (p$class [svms] * (p$data [svms,] %*% w))^2;
      circle.z <- seq (-1, 1, length=10000);
      bound <- - w [1]^2 * (4 * w [2]^2 * circle.z^2 -
                            4 * w [2]^2 +
                            4 +
                            4 * w [3]^2 * circle.z^2 -
                            8 * w [3] * circle.z +
                            4 * w [3] * circle.z * gamma -
                            4 * gamma +
                            gamma * gamma +
                            4 * circle.z^2 * w [1]^2 -
                            4 * w [1]^2);

      circle.x       <- rep (-10, length (circle.z));
      circle.y       <- rep (-10, length (circle.z));
      idx            <- (bound >= 0);
      circle.y [idx] <- -0.5 * (-2 * w [2] +
                                2 * w [2] * w [3] * circle.z [idx] +
                                gamma * w [2] -
                                sqrt (bound [idx])) / (w [1]^2 + w [2]^2);
      circle.x [idx] <- -0.5 * (w [2] * sqrt (bound [idx]) +
                                2 * w [1]^2 * w [3] * circle.z [idx] -
                                2 * w [1]^2 +
                                w [1]^2 * gamma) / (w [1]^2 + w [2]^2) / w [1];
      lines (project.3d (circle.x [idx], circle.y [idx], circle.z [idx],
                         VT=VT),
             col="black", lwd=4);

      circle.y [idx] <- -0.5 * (-2 * w [2] +
                                2 * w [2] * w [3] * circle.z [idx] +
                                gamma * w [2] +
                                sqrt (bound [idx])) / (w [1]^2 + w [2]^2);
      circle.x [idx] <- -0.5 * (-w [2] * sqrt (bound [idx]) +
                                2 * w [1]^2 * w [3] * circle.z [idx] -
                                2 * w [1]^2 +
                                w [1]^2 * gamma) / (w [1]^2 + w [2]^2) / w [1];
      lines (project.3d (circle.x [idx], circle.y [idx], circle.z [idx],
                         VT=VT),
             col="black", lwd=4);
    }
  }
  
  return (p);
}


######################################################################
##  F I S H E R   L I N E A R   D I S C R I M I N A N T S
######################################################################

## classification with Fisher linear discriminants
perc.fisher <- function (p,             # old linear classifier
                         lambda=1e-3,   # regularisation variance
                         verbosity=0) { # output convergence msgs?
  
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  if (p$classifier == FALSE)
    stop ("perc.GP.class can only be applied to regression");
  if (lambda < 0)
    stop ("Lambda has to be positive");

  ## extract training set size
  m <- length (p$alpha);

  ## compute the cardinality of the two class sets
  m.plus <- sum (p$class == +1);
  m.minus <- m - m.plus;

  ## construct the temporary g-vectors
  g.plus <- cbind (rep (0, length=m));
  g.plus [p$class == +1] <- 1/m.plus;
  g.plus <- p$K %*% g.plus;
  g.minus <- cbind (rep (0, length=m));
  g.minus [p$class == -1] <- 1/m.minus;
  g.minus <- p$K %*% g.minus;

  ## compute the matrix to be inverted
  delta <- g.plus - g.minus;
  H <- 1/m * (p$K %*% p$K - m.plus * g.plus %*% t (g.plus) -
              m.minus * g.minus %*% t (g.minus)) +
                diag (rep (lambda, length=m));
  p$alpha <- solve (H, delta);
  
  ## create the vector and diagonal matrix of ones
  ## type conversion only
  p$alpha <- as.vector (p$alpha);
  p$margin <- 0;

  H.inv <- solve (H);
  p$theta <- as.vector (1/2 * (t (g.minus) %*% H.inv %*% g.minus -
                               t (g.plus) %*% H.inv %*% g.plus) +
                        log (m.plus / m.minus));
  
  return (p);
}

############################################################
##     P A R Z E N   W I N D O W S
############################################################

perc.parzen <- function (p) {
  if (class (p) != "perc")
    stop ("p (1st arg) has to be a perceptron");
  
  p$alpha <- p$class;        # set all alpha's to +1/-1 (for RBF kernels
  p$theta <- 0;              # this is equivalent to kernel density estimation
  p$margin <- 0;             # of the class conditional densities
  return(p);
}

############################################################
##     P L O T T I N G   R O U T I N E S
############################################################


############################################################
## computes an appropriate palette
############################################################

perc.colors <- function (n) {
  k <- n / 2;
  c (hsv (h = 9/12, s = seq (0.5, 0, length = k + 1) [-k - 1], v = 1),
     hsv (h = 2/12, s = seq (0, 0.5, length = n - k + 1) [-1], v = 1));
}

############################################################
## plot the class regions (in 2D) and exit in higher dimensions
############################################################

plot.perc <- function (p, nocol=128, grid=50,
		       mag=.2, image=TRUE, margin=TRUE, mark=TRUE,
                       pts=TRUE, dec.col="black", mar.col="black",
                       ...) {

  ##
  ## this block is for plotting classifiers
  ## 

  if (p$classifier == TRUE) {
    if (dim (p$data) [2] == 2 ||
        (dim (p$data) [2] == 3 && (min (p$data [,3]) == max (p$data [,3])))) {
      olwd <- par()$lwd;
      par (lwd=2);
      xl     <- range (p$data [,1]);
      xlen   <- xl [2] - xl [1];
      xl [1] <- xl [1] - mag*xlen;
      xl [2] <- xl [2] + mag*xlen;
      yl     <- range (p$data [,2]); 
      ylen   <- yl [2] -yl [1];
      yl [1] <- yl [1] - mag*ylen;
      yl [2] <- yl [2] + mag*ylen;
      x      <- seq (xl [1], xl [2], length=grid);
      y      <- seq (yl [1], yl [2], length=grid);
      Y      <- rep (y, rep (length (x), length (y)));
      X      <- rep (x, length.out = length (Y));
      pt     <- cbind (X,Y);
      if (dim (p$data) [2] == 3) {
        z    <- perc.classify (cbind (pt, p$data [1,3]), p);
      } else {
        z    <- perc.classify (pt, p);
      }
      
      ## print out an image 
      if (image) {
        image (x, y, matrix (z, length (x)),
               col=perc.colors (nocol),
               xlab="X",
               ylab="Y", ...);
        contour (x, y, matrix (z, length (x)),
                 level=0, col=dec.col, add=T);
        if (p$margin > 0 && margin) {
          contour (x, y, matrix (z, length (x)),
                   level=p$margin, lty=2, col=mar.col, add=T);
          contour (x, y, matrix (z, length (x)),
                   level=-p$margin, lty=2, col=mar.col, add=T);
        }
        if (pts == TRUE) {
          points (p$data [p$class==1, 1], p$data[p$class==1, 2],
                  col="red", pch = 19);
          points (p$data [p$class==-1,1], p$data[p$class==-1,2],
                  col="blue", pch = 19);
          if (p$eps > 0 && mark) {
            points (p$data [abs (p$alpha) > p$eps, 1],
                    p$data [abs (p$alpha) > p$eps, 2], col="black",
                    pch = 1, cex=1.5);
          }
        }
      }
      ## print out an b&w plot
      else {
        contour (x, y, matrix (z, length (x)),
                 level=0, col=dec.col, xlab="X", ylab="Y", ...);
        if (p$margin > 0 && margin) {
          contour (x, y, matrix (z, length (x)),
                   level=p$margin, lty=2, col=mar.col, add=T);
          contour (x, y, matrix (z, length (x)),
                   level=-p$margin, lty=2, col=mar.col, add=T);
        }
        if (pts == TRUE) {
          points (p$data [p$class==1,1], p$data [p$class==1,2],
                  col="black", pch = 20,cex=1);
          points (p$data [p$class==-1,1], p$data [p$class==-1,2],
                  col="black", pch = 4, cex=1);
          if (p$eps > 0 && mark) {
            points (p$data [abs(p$alpha)>p$eps,1],
                    p$data [abs(p$alpha)>p$eps,2],
                    col="black", pch = 1, cex=2);
          }
        }
      }
      par (lwd = olwd);
    }
    else {
      stop("Plots are only for 2D perceptrons available.");
    }
  } else {
    
    ##
    ## this block is for plotting regression functions
    ## 
    
    if (dim (p$data) [2] == 1) {
      olwd <- par()$lwd;
      par (lwd=2);
      xl     <- range (p$data [,1]);
      xlen   <- xl [2] - xl [1];
      xl [1] <- xl [1] - mag*xlen;
      xl [2] <- xl [2] + mag*xlen;
      x      <- seq (xl [1], xl [2], length=grid);
      y      <- perc.classify (cbind (x), p);
      yl     <- range (c (y-p$margin, y+p$margin, p$class));
      ylen   <- yl [2] - yl [1];
      yl [1] <- yl [1] - mag*ylen;
      yl [2] <- yl [2] + mag*ylen;
      
      ## print out an image 
      plot (x, y, type="l", xlim=xl, ylim=yl,
            col="black", xlab="X", ylab="Y", ...);
      if (margin) {
        lines (x, y+p$margin, type="l", lty=2, col="black");
        lines (x, y-p$margin, type="l", lty=2, col="black");
      }
      if (pts == TRUE) {
        points (p$data, p$class, col="red", pch = 19);
        if (mark) {
          points (p$data [abs (p$alpha) > p$eps],
                  p$class [abs (p$alpha) > p$eps],
                  col="black", pch = 1, cex=1.5);
        }
      }
      par (lwd = olwd);
    } else {
      stop("Plots are only for 1D regression perceptrons available.");
    }
  }
}

############################################################
## this set of functions if for 3D plotting (for version space
## pictures)
############################################################

plot.3d <- function (x, y, z,        # the list of 3D coordinates 
                     theta=30,       # the x-axis angle
                     phi=60,         # the y-axis angle
                     scale=1,        # the scale of the eye distance
                     radius,         # either scale of radius is used (see above)
                     pch=1,          # charachter to be used for plots
                     lwd=1,          # line width 
                     cex=1,          # character expansion
                     col="black",    # color for plots
                     ...) { 
  if (missing (z) && missing (y)) {
    if (!is.matrix (x)) {
      stop ("1st argument has to be a matrix");
    } else {
      X <- x;
    }
  } else {
    X <- cbind (x, y, z);
  }

  ##  if (theta < 0 || theta > 360)
  ##    stop ("theta has to be in the range [0:360]");
  ##  if (phi < 0 || phi > 360)
  ##    stop ("phi has to be in the range [0:360]");
  if (!missing (radius)) {
    R <- radius;
  } else {
    R <- 2 * max (abs (X)) * scale;
  }
  theta  <- theta * pi / 180;
  phi    <- phi * pi / 180;
  
  Xrot <- rbind (c (1        , 0           , 0          ),
                 c (0        , cos (theta) , sin (theta)),
                 c (0        , -sin (theta), cos (theta)));
  Yrot <- rbind (c (cos (phi), 0           , -sin (phi) ),
                 c (0        , 1           , 0          ),
                 c (sin (phi), 0           , cos (phi)  ));
  XYrot <- Xrot %*% Yrot;
  
  ## compute the projection
  VT <- list (rot=XYrot, trans=-R);
  Xp <- project.3d (X, VT=VT);
  plot (range (Xp [,1]), range (Xp [,2]), type = "n", ...);
  points (Xp [,1], Xp [,2], pch=pch, lwd=lwd, col=col, cex=cex);
  
  return (VT);
}

############################################################
## projects a given set of 3D points onto a plane
############################################################

project.3d <- function (x, y, z,        # the list of 3D coordinates 
                        VT,             # the return of plot.3d
                        hidden=FALSE) { # clears hidden points 
  if (missing (z) && missing (y)) {
    if (!is.matrix (x)) {
      stop ("1st argument has to be a matrix");
    } else {
      X <- x;
    }
  } else {
    X <- cbind (x, y, z);
  }

  Xtmp <- X %*% VT$rot;
  Xp <- cbind (Xtmp [,1],
               Xtmp [,2]);
  
  if (hidden) {
    return (Xp [Xtmp [,3] < 0,]);
  } else {
    return (Xp);
  }
}


############################################################
## computes the grid points of a 3D ball
############################################################

ball <- function (radius=1,       # radius of the ball 
                  ntheta=36,      # number of steps in x-angle 
                  nphi=36) {      # number of steps in y-angle
  X <- matrix (0, nrow=ntheta*nphi, ncol=3);
  i <- 1;
  for (theta in seq (0, 2*pi, length=ntheta)) {
    for (phi in seq (pi/2, 3*pi/2, length=nphi)) {
      X [i,] <- radius * c (-sin (phi),
                            sin  (theta) * cos (phi),
                            cos  (theta) * cos (phi));
      i <- i + 1;
    }
  }
  return (X);
}

############################################################
## this routine extents the standard errorbar plot routines
############################################################

plot.errorbars <- function (x, y, error,
                            lwd=1,        # width of main curve
                            lty=1,        # style of main curve
                            col="black",  # color of main curve
                            blwd=lwd,     # line width of bars
                            blty=lty,     # line style of bars
                            bcol=col,     # color of bars
                            ticks=TRUE,   # draw ticks on the errorbars?
                            ticch=-1,     # the tic-character (-1 means bars)
                            add=FALSE,    # should the plot be added to a
                                          # existing one?
                            ...) { 
  if (!is.vector (x))
    stop ("1st arg (x) has to be a vector");
  if (!is.vector (y))
    stop ("2nd arg (y) has to be a vector");
  if (!is.vector (error))
    stop ("3rd arg (y) has to be a vector");
  if (length (x) != length (y) || length (x) != length (error)) 
    stop ("all arguments should have the same length");

  ## determine the range for outputs
  yl <- range (c (y+error, y-error));
  yd <- yl [2] - yl [1]; yl [1] <- yl [1] - .05 * yd;
  yl [2] <- yl [2] + .05 * yd;
  xl <- range (x);
  ## compute optimal tick length
  tl <- (xl [2] - xl [1])/(10 * length (x));
  
  ## plot a template which fits at least the errorbars
  if (!add) {
    plot (xl, yl, type="n", ...);
  }
  ## plot the lines through (x,y)
  lines (x, y, type="l", lwd=lwd, lty=lty, col=col);
  
  ## main loop for the error bars
  for (i in 1:length (error)) {
    lines (c (x [i], x [i]), c (y [i] - error [i], y [i] + error [i]),
           lwd=blwd, lty=blty, col=bcol);
    if (ticks) {
      if (ticch == -1) {
        lines (c (x [i] - tl, x [i] + tl), c (y [i] - error [i],
                                              y [i] - error [i]),
               lwd=blwd, lty=blty, col=bcol);
        lines (c (x [i] - tl, x [i] + tl), c (y [i] + error [i],
                                              y [i] + error [i]),
               lwd=blwd, lty=blty, col=bcol);
      } else {
        points (c (x [i], x [i]), c (y [i] - error [i], y [i] + error [i]),
                pch=ticch, lwd=blwd, col=bcol);
      }
    }
  }
}

############################################################
##     D A T A S E T   R O U T I N E S
############################################################

############################################################
## reads a data set from a file
############################################################

read.dataset <- function (file="data", verbosity=0) {
  if (!file.exists (file))
    stop ("File does not exists!");
  header <- scan (file, n = 3, quiet=!verbosity);
  if (header [3] != 0)
    stop ("File is of incorrect format!");
  d <- header [1];
  n <- header [2];
  D <- matrix (scan (file, skip=1, n=(d+1)*n, quiet=!verbosity),
               nrow=n, byrow=T);
  return (list (data=D[,1:d], class=D[,(d+1)]));
}

############################################################
## normalise the dataset in input space
############################################################

normalize.dataset <- function (X, verbosity=0) {
  if (!is.matrix (X))
    stop ("X (1st arg) has to be a matrix");

  if (verbosity == 2)
    cat ("[normalize dataset]\n");
  for (i in 1:dim(X) [1]) {
    X [i,] <- X [i,] / sqrt (X[i,] %*% X [i,]);
  }
  return (X);
}

############################################################
##     E X A M P L E   C O D E
############################################################

exmp.bayes <- function () {
  load ("linear/artifical4.dat");
  p <- perc.new (cbind (data$x, data$y, 1), data$class, kernel.new ("linear"));
  p <- perc.svm (p);
  ## p <- perc.learn (p);
  pb <- perc.bayes (p, max.epoch = 1000, verb=2, plot=T, 
                    theta=60, phi=-30, radius=1.5);
  
  par (mfrow = c(1, 2));
  plot (p, margin=F, mark=F, nocol=32); title ("SVM");
  plot (pb, margin=F, mark=F, nocol=32); title ("BPM");
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.bayes2 <- function () {
  load ("linear/artifical1.dat");
  p <- perc.new (cbind (data$x, data$y), data$class, kernel.new ("rbf", sigma=1));
  p <- perc.svm (p);
  pb <- perc.bayes (p, disp.epoch=500, max.epoch = 1000, verb=2);
  
  plot (pb, margin=F, mark=F); title ("BPM");
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.svm <- function (file="linear/artifical1.dat", C=1e20,
                      kernel=kernel.new ("rbf", sigma=1), norm=FALSE) {
  load (file);
  p <- perc.new (cbind (data$x, data$y), data$class, kernel);
  p <- perc.svm (p, C=C, norm=norm);
  
  par (mfrow = c(1, 1));
  plot (p); title ("SVM");
  cat ("margin : ", formatC (perc.new.margin (p)), "\n")
  readline ("Press <Return> to finish.");
}

exmp.GP.class <- function (file="linear/artifical1.dat", C=1,
                           kernel=kernel.new ("rbf", sigma=1),
                           tol=1e-4, var=0.1) {
  load (file);
  p <- perc.new (cbind (data$x, data$y), data$class, kernel);
  p <- perc.GP.class (p, C=C, var=var, tol=tol);
  
  par (mfrow = c(1, 1));
  plot (p); title ("GPC");
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.GP.regress <- function (file="linear/regression1.dat", var=0.1,
                             kernel=kernel.new ("rbf", sigma=1)) {
  load (file);
  p <- perc.new (cbind (data$x), data$y, kernel);
  p <- perc.GP.regress (p, var=var);
  
  par (mfrow = c(1, 1));
  plot (p, grid=1000); title ("GPR");
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.RVM.regress <- function (file="linear/regression1.dat", tol=.Machine$double.eps,
                              kernel=kernel.new ("rbf", sigma=1), ...) {
  load (file);
  p <- perc.new (cbind (data$x), data$y, kernel);
  p <- perc.RVM.regress (p, tol=tol, verbosity=1, iterations=100, ...);
  
  par (mfrow = c(1, 1));
  plot (p); title ("RVM");
  readline ("Press <Return> to finish.");
  dev.off ()
}

exmp.RVM.class <- function (file="linear/artifical1.dat",
                            kernel=kernel.new ("rbf", sigma=1),
                            C=1,
                            tol=1e-4, verbosity=2, ...) {
  load (file);
  p <- perc.new (cbind (data$x, data$y), data$class, kernel);
  p <- perc.RVM.class (p, C=C, ...);
  
  par (mfrow = c(1, 1));
  plot (p); title ("RVMC");
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.fisher <- function (file="linear/artifical1.dat", lambda=1e-3,
                               kernel=kernel.new ("rbf", sigma=1), ...) {
  load (file);
  p <- perc.new (cbind (data$x, data$y), data$class, kernel);
  p <- perc.fisher (p, lambda=lambda, ...);
  
  par (mfrow = c(1, 1));
  plot (p); title ("KFD");
  cat ("Classification error:", perc.class.error (p$data, p$class, p));
  readline ("Press <Return> to finish.");
  dev.off ();
}

exmp.perc <- function (margin=0, method="standard",
                       kernel=kernel.new ("rbf", sigma=1)) {
  load ("linear/artifical5.dat");
  p <- perc.new (cbind (data$x, data$y), data$class, kernel);
  p <- perc.learn (p, margin=margin, method=method, verbosity=1);
  
  par (mfrow = c(1, 1));
  plot (p); title (paste ("Dual Perceptron (Margin=", formatC (margin), ")",
                          sep=""));
  cat ("margin : ", formatC (perc.new.margin (p)), "\n")
  readline ("Press <Return> to finish.");
}

exmp.nu.svm <- function () {
  load ("linear/artifical1.dat");
  p <- perc.new (cbind (data$x, data$y), data$class, kernel.new ("rbf", sigma=1));
  p2 <- perc.nu.svm (p, nu=.2);

  par (mfrow = c(1, 1));
  plot (p2);
  title (expression (paste (nu, "-SVM")));
  readline ("Press <Return> for next plot.");
  
  par (mfrow=c (3,3))
  for (i in c (.01, .02, .03, .04, .05, .06, .07, .08, .09)) {
    p2 <- perc.nu.svm (p, nu=i);
    plot(p2);
    title (expression (paste (nu, "=", i)));
  }
  readline ("Press <Return> to finish.");
  dev.off ();
}
  
exmp.plot.3d <- function (n = 50) {
  ##  vs <- matrix (scan ("ball.dat"), ncol=3, byrow=TRUE);
  b <- ball (ntheta=n, nphi=n);
  for (theta in seq (0, 180, by=20)) {
    for (phi in seq (0, 360, by=20)) {
      VT <- plot.3d (b, theta=theta, radius=10, phi=phi, pch='.');
      for (i in 1:n) {
        p <- project.3d (b [seq (i, n*n, by = n),], VT=VT);
        lines (p [,1], p [,2], col="red", pch=20);
      }
      for (i in 1:n) {
        p <- project.3d (b [seq ((i-1)*n+1,i*n),], VT=VT);
        lines (p [,1], p [,2], col="blue", pch=20);
      }
      ##      vs.proj <- project.3d (vs, VT=VT,hidden=T);
      ##      points (vs.proj, col="black", pch="*");
      cat ("theta: ", theta, "\tphi: ", phi, "\n");
      readline ("Hit <Return> for next plot");
    }
  }
}
