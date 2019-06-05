###
### demo        generates plots of several distributions
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
## density of the discrete uniform distribution 
############################################################

duniform <- function (i, min=0, max=1) {
  ret <- rep (0, length (i));
  ret [i >= min & i <= max] <- 1/(max - min + 1);
  
  return (ret);
}

############################################################
## DISCRETE probability distributions
############################################################

plot.discrete <- function (output) {
  ## Bernoulli distribution
  i <- 0:1;

  X <- rbind (dbinom (i, size=1, prob=0.2),
              dbinom (i, size=1, prob=0.5),
              dbinom (i, size=1, prob=0.9));
  rownames (X) <- c ("p=0.2", "p=0.5", "p=0.9");
  colnames (X) <- c ("0", "1");

  if (output == 'PS') {
    postscript (file="distribution_bernoulli.ps", horizontal=TRUE);
  }
  par (cex.axis=1.5, cex.lab=2.0, cex=1.4);
  par (mai=c (1.25,1.35,0.25,0.25));
  barplot (X, beside=TRUE,
           legend.text=rownames (X),
           col = c (gray (0.6), gray (0.4), gray (0.9)),
           ylab = expression (paste (bold (P), group ("(", "X=i", ")"))),
           xlab = "i");
  if (output == 'PS') {
    cat ("[Bernoulli plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
  
  ## Binomial distribution
  i <- 0:20;

  X <- rbind (dbinom (i, size=20, prob=0.2),
              dbinom (i, size=20, prob=0.5),
              dbinom (i, size=20, prob=0.9));
  rownames (X) <- c ("p=0.2", "p=0.5", "p=0.9");
  colnames (X) <- c ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19",
                     "20");

  if (output == 'PS') {
    postscript (file="distribution_binomial.ps", horizontal=TRUE);
  }
  par (cex.axis=1.5, cex.lab=2.0, cex=1.4);
  par (mai=c (1.25,1.35,0.25,0.25));
  barplot (X, beside=TRUE,
           legend.text=rownames (X),
           col = c (gray (0.6), gray (0.4), gray (0.9)),
           ylab = expression (paste (bold (P), group ("(", "X=i", ")"))),
           xlab = "i");
  if (output == 'PS') {
    cat ("[Binomial plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();

  ## Poisson distribution
  i <- 0:20;

  X <- rbind (dpois (i, lambda=1),
              dpois (i, lambda=2),
              dpois (i, lambda=10));
  leg <- c (expression (paste (lambda, "=1")),
            expression (paste (lambda, "=2")),
            expression (paste (lambda, "=10")));
  colnames (X) <- c ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19",
                     "20");

  if (output == 'PS') {
    postscript (file="distribution_poisson.ps", horizontal=TRUE);
  }
  par (cex.axis=1.5, cex.lab=2.0, cex=1.4);
  par (mai=c (1.25,1.35,0.25,0.25));
  barplot (X, beside=TRUE,
           legend.text = leg,
           col = c (gray (0.6), gray (0.4), gray (0.9)),
           ylab = expression (paste (bold (P), group ("(", "X=i", ")"))),
           xlab = "i");
  if (output == 'PS') {
    cat ("[Poisson plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();

  ## Uniform distribution
  i <- 0:20;

  X <- rbind (duniform (i, min=0, max=1),
              duniform (i, min=0, max=5),
              duniform (i, min=0, max=20));
              
  rownames (X) <- c ("b=1", "b=5", "b=20");
  colnames (X) <- c ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19",
                     "20");

  if (output == 'PS') {
    postscript (file="distribution_uniform.ps", horizontal=TRUE);
  }
  par (cex.axis=1.5, cex.lab=2.0, cex=1.4);
  par (mai=c (1.25,1.35,0.25,0.25));
  barplot (X, beside=TRUE,
           legend.text = rownames (X),
           col = c (gray (0.6), gray (0.4), gray (0.9)),
           ylab = expression (paste (bold (P), group ("(", "X=i", ")"))),
           xlab = "i");
  if (output == 'PS') {
    cat ("[Uniform plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
}

############################################################
## CONTINOUS probability distributions
############################################################

plot.continous <- function (output) {
  ## Normal distribution
  x <- seq (-5, 5, length=200);

  X <- rbind (dnorm (x, mean=0, sd=1),
              dnorm (x, mean=1, sd=0.5),
              dnorm (x, mean=0, sd=2));
  
  if (output == "PS") {
    postscript (file="distribution_normal.ps", horizontal=TRUE);
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (c (-5, 5), range (X), type="n",
        ylab = expression (paste (bold (f), group ("(", "x", ")"))),
        xlab="x", cex.lab=3.5, cex.axis=2);
  lines (x, X [1,], lwd=4, col=gray (0.1));
  lines (x, X [2,], lwd=4, col=gray (0.4));
  lines (x, X [3,], lwd=4, col=gray (0.7));

  text (x [90],  X [1,95],
        expression (paste (mu, "=0, ", sigma, "=1")),   cex=2.5, pos=3);
  text (x [122], X [2,122], 
        expression (paste (mu, "=1, ", sigma, "=0.5")), cex=2.5, pos=4);
  text (x [50],  X [3,50],
        expression (paste (mu, "=0, ", sigma, "=2")),   cex=2.5, pos=2);
  if (output == 'PS') {
    cat ("[Normal plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
  
  ## Exponential distribution
  x <- seq (0, 5, length=200);

  X <- rbind (dexp (x, rate=1),
              dexp (x, rate=2),
              dexp (x, rate=10));

  if (output == 'PS') {
    postscript (file="distribution_exponential.ps", horizontal=TRUE);
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (c (0, 5), c (0, 1), type="n",
        ylab = expression (paste (bold (f), group ("(", "x", ")"))),
        xlab="x", cex.lab=3.5, cex.axis=2);
  lines (x, X [1,], lwd=4, col=gray (0.0));
  lines (x, X [2,], lwd=4, col=gray (0.4));
  lines (x, X [3,], lwd=4, col=gray (0.7));

  text (x [100], X [1,100]+0.05,
        expression (paste (lambda, "=1")),  cex=2.5, pos=4);
  text (x [23], X [2,23], 
        expression (paste (lambda, "=2")),  cex=2.5, pos=4);
  text (x [22], X [3,22],
        expression (paste (lambda, "=10")), cex=2.5, pos=4);
  if (output == 'PS') {
    cat ("[Exponential plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
  
  ## Gamma distribution
  x <- seq (0, 1.2, length=200);

  X <- rbind (dgamma (x, shape=2, scale=0.1),
              dgamma (x, shape=5, scale=0.1),
              dgamma (x, shape=1, scale=0.25));

  if (output == 'PS') {
    postscript (file="distribution_gamma.ps", horizontal=TRUE);
  }  
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (range (x), c (0, 4), type="n",
        ylab = expression (paste (bold (f), group ("(", "x", ")"))),
        xlab="x", cex.lab=3.5, cex.axis=2);
  lines (x, X [1,], lwd=4, col=gray (0.0));
  lines (x, X [2,], lwd=4, col=gray (0.4));
  lines (x, X [3,], lwd=4, col=gray (0.7));

  text (x [30], X [1,30],
        expression (paste (alpha, "=2, ", beta, "=0.1")),  cex=2.5, pos=4);
  text (x [70], X [2,70] + 0.05, 
        expression (paste (alpha, "=5, ", beta, "=0.1")),  cex=2.5, pos=4);
  text (x [90], X [3,90] + 0.05,
        expression (paste (alpha, "=1, ", beta, "=0.25")),  cex=2.5, pos=4);
  if (output == 'PS') {
    cat ("[Gamma plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
  
  ## Beta distribution
  x <- seq (0, 1, length=200);

  X <- rbind (dbeta (x, shape1=0.5, shape2=0.5),
              dbeta (x, shape1=1, shape2=1),
              dbeta (x, shape1=5, shape2=5),
              dbeta (x, shape1=2, shape2=5));

  if (output == 'PS') {
    postscript (file="distribution_beta.ps", horizontal=TRUE);
  }  
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (range (x), c (0, 2.6), type="n",
        ylab = expression (paste (bold (f), group ("(", "x", ")"))),
        xlab="x", cex.lab=3.5, cex.axis=2);
  lines (x, X [1,], lwd=4, col=gray (0.0));
  lines (x, X [2,], lwd=4, col=gray (0.3));
  lines (x, X [3,], lwd=4, col=gray (0.6));
  lines (x, X [4,], lwd=4, col=gray (0.8));

  text (x [85], X [1,85] - 0.05,
        expression (paste (alpha, "=0.5, ", beta, "=0.5")),  cex=2.5, pos=1);
  text (x [120], X [2,120], 
        expression (paste (alpha, "=1, ", beta, "=1")),  cex=2.5, pos=3);
  text (x [100], X [3,100],
        expression (paste (alpha, "=5, ", beta, "=5")),  cex=2.5, pos=3);
  text (x [40], X [4,40],
        expression (paste (alpha, "=2, ", beta, "=5")),  cex=2.5, pos=3);
  if (output == 'PS') {
    cat ("[Beta plot generated]\n");
  } else {
    readline("Press any key for next plot");
  }
  dev.off ();
  
}

############################################################
## NORMAL DENSITY plot
############################################################

plot.normal <- function (output) {
  ## Normal distribution
  x <- seq (-3, 3, length=50);
  y <- seq (-3, 3, length=50);

  ## define the orientation and excentritiy of the covariance ellipses
  sigma.x <- c (1, 1/2);
  sigma.y <- c (1, 1);
  rho <- c (0, 0.5);

  for (j in 1:length (rho)) {
    ## compute the density 
    z <- outer (x, y, function (x, y) {
      
      ## compute the covariance matrix
      C <- diag (c (sigma.x [j], sigma.y [j])) %*%
        cbind (c (1, rho [j]), c (rho [j], 1)) %*%
          diag (c (sigma.x [j], sigma.y [j]));

      ## invert this matrix and compute the normalisation constant
      norm <- 1/(2 * pi * sqrt (det (C)));
      ret <- double (0);
      C.inv <- solve (C);
      
      ## compute the plot point by point
      for (i in 1:length (x)) {
        ret <- c (ret, norm * exp (-1/2 *
                                   rbind (c (x [i], y [i])) %*% C.inv %*%
                                   cbind (c (x [i], y [i]))));
      }
      return (ret);
    });
  
    if (output == 'PS') {
      postscript (file=paste ("multi_normal", j, ".ps", sep=""),
                  horizontal=TRUE);
    }
    persp (x, y, z, theta=50, phi=20, shade=0.2,
           lwd=0.5, ticktype="detailed", xlab="x", ylab="y",
           zlab="density", cex.lab=1.2);
    
    if (output == 'PS') {
      cat ("[Multidimensional normal no", j, "plot created]\n");
    } else {
      readline("Press any key for next plot");
    }
    dev.off ();
  }
}

############################################################
## BOOK code
############################################################

book <- function (output='SCREEN') {
  plot.discrete (output);
  plot.continous (output);
  plot.normal (output);
}

book()
