### demo        creates all plots for the covering number section
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## little dataset (each row is a datum)
############################################################

X <- c (1, 2);

############################################################
## generate the cloud of points that needs to be
## covered using the RBF kernel
############################################################

gen.points <- function (X, N = 100, sigma=0.5) {
  
  ## compute the inner product between x1 and x2
  k <- exp (-(X [1] - X [2])^2/(2*sigma^2));

  ## generate Gram matrix
  G <- cbind (c (1, k), c (k, 1));
  
  ## generate admissible alpha's
  a <- runif (N, min=-sqrt (1/(1-k^2)), max=sqrt (1/(1-k^2)));
  alpha <- cbind (c (a), runif (N, min=(-a * k - sqrt (a^2 * (k^2 - 1) + 1)),
                                max=-a * k + sqrt (a^2 * (k^2 - 1) + 1)));
  
  out <- t (G %*% t (alpha));
  return (list (alpha=alpha, out=out));
}

############################################################
## computes the function values for a specific 
## dataset and RBF kernel
############################################################

draw.fct <- function (X, alpha, x, sigma=0.5) {
  
  ## compute the inner product between x1 and x2
  k1 <- exp (-(X [1] - x)^2/(2*sigma^2));
  k2 <- exp (-(X [2] - x)^2/(2*sigma^2));
  H <- cbind (k1, k2);
  y <- H %*% alpha;

  return (y);
}

############################################################
## compute a cover of the function class (greedily)
############################################################

cover <- function (out, gamma=0.5) {

  cov <- double (0);
  
  while (nrow (out) > 0) {

    ## compute the coverage of all points
    coverage <- double (0);
    for (i in 1:nrow (out)) {
      ## remove all points that are covered
      rem.idx.x <- (abs (out [, 1] - out [i, 1]) < gamma);
      rem.idx.y <- (abs (out [, 2] - out [i, 2]) < gamma);
      coverage <- c (coverage, sum (rem.idx.x & rem.idx.y));
    }

    ## choose the point with the highest coverage
    idx <- order (coverage) [1];

    ## remove all covered points
    cov <- rbind (cov, out [idx, ]);
    rem.idx.x <- (abs (out [, 1] - out [idx, 1]) < gamma);
    rem.idx.y <- (abs (out [, 2] - out [idx, 2]) < gamma);
    out <- out [!(rem.idx.x & rem.idx.y), ];
  }

  return (cov);
}



############################################################
## BOOK plots (set output to 'PS' for PostScript files)
############################################################

book <- function (output='SCREEN') {

  ## initialise "pseudo" random number generator
  set.seed (42);
  
  ## generate a dataset and cover
  sigma <- 0.3;
  gamma <- 0.1;
  R <- gen.points (X, N=200, sigma=sigma);
  cov <- cover (R$out, gamma=gamma);

  ############################################################
  ## draw a few functions
  if (output == 'NEW_PS') {
    if (!file.exists ("cover_fct1.ps")) {
      postscript (file="cover_fct1.ps");
    }
  }
  
  par (mai=c (1.0,1.25,0.25,0.25));
  N <- 20;
  d <- max (X) - min (X);
  x <- seq (min (X) - d * 0.1, max (X) + d * 0.1, length=100);
  y <- double (0);
  ylim <- c (0, 0);

  for (i in 1:N) {
    tmp <- matrix (draw.fct (X, R$alpha [i, ], x, sigma=sigma), nrow=1);
    ylim <- range (ylim, tmp);
    y <- rbind (y, tmp);
  }
  
  ## draw little dataset with some functions
  plot (range (x), ylim, type="n",
        cex.lab=3.5, cex.axis=2.5,
        xlab="x", ylab="f(x)");
  for (i in 1:N) {
    lines (x, y [i, ], lwd=2);
  }
  lines (x, rep (0, length=length (x)), lwd=1);
  points (X, rep (0, length=length (X)), pch=4, lwd=4, cex=3);
      
  readline ("Press any key for next plot");
  ## dev.off ();

  ############################################################
  ## draw a cover
  if (output == 'NEW_PS') {
    if (!file.exists ("cover_1.ps")) {
      postscript (file="cover_1.ps");
    }
  }
    
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (c (R$out [, 1] + gamma, R$out [, 1] - gamma),
        c (R$out [, 2] + gamma, R$out [, 2] - gamma),
        type="n", cex.lab=3.5, cex.axis=2.5,
        xlab=expression (paste ("            ", f(x[1]))),
        ylab=expression (f(x[2])));
  for (i in 1:nrow (cov)) {
    rect (cov [i, 1] - gamma, cov [i, 2] - gamma,
          cov [i, 1] + gamma, cov [i, 2] + gamma, lwd=2, col=gray(0.9));
    points (cov [i, 1], cov [i, 2], pch=19, cex=2);
  }
  for (i in 1:nrow (R$out)) {
    if (sum (R$out [i, 1] == cov [, 1]) +
        sum (R$out [i, 2] == cov [, 2]) == 0) {
      points (R$out [i, 1], R$out [i, 2], pch=19, col=gray(0.4));
    }
  }
  
  readline ("Press any key for next plot");
  ## dev.off ();

}

book ();


