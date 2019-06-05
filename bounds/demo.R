###
### demo.R           1) depicts the growth of x*(1-ln(x))
###                  2) depicts the growth of the fat bound
###                  3) depicts the growth of the margin bound
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
## computes the valued of the fat bound (ignoring delta)
############################################################
fat.bound <- function (m, d=1, error=1) {
  return (2/m * (1+ceiling (d * log(8*exp(1)*m/d,2) * log(32*m,2))) -
          error);
}

############################################################
## computes the derivative of the fat bound (for fixed d)
############################################################
fat.bound.derivative <- function (m, d=1) {
  return (-2 * (1+(d*log(8*exp(1)*m/d)*log(32*m)/log(2)^2))/m^2
          +2 * (d*log(32*m)/(m*log(2)^2) +
                d*log(8*exp(1)*m/d)/(log(2)^2*m))/m);
}

############################################################
## computes the minimal value of m for a given d
## in the fat bound case
############################################################
fat.min.training.set <- function (d, error=1, verbosity=0) {
  min.m <- integer (0);
  for (i in 1:length (d)) {
    m <- 100 * d [i];
    if (verbosity > 0) {
      cat ("Determining the minimal training set size for d=", d [i], "\n");
    }
    while (abs (fat.bound (m, d [i], error)) > 1e-4) {
      if (verbosity > 0) {
        cat ("m=", m, "\tbound=", fat.bound (m, d [i], error), "\n");
      }
      m <- m - fat.bound (m, d [i], error)/fat.bound.derivative (m, d [i]);
    }
    if (verbosity > 0) {
      cat ("convergence suceeded\n\n");
    }
    min.m <- c (min.m, floor (m));
  }
  return (min.m);
}

############################################################
## plots the minimal required training set size
############################################################

fat.make.plot <- function (d, error=1) {
  
  m <- fat.min.training.set (d, error=error);

  plot (d, m, type="n",
        xlab=expression (paste ("effective complexity = ",
            fat, bgroup("(", frac(gamma,8), ")"))),
            ylab="minimal training set size",
        cex.axis = 2.5, cex.lab=3.5);
  lines (d, m, lwd=4);
  points (d, m, pch=19, lwd=2, cex=2);
}

############################################################
## computes the valued of the bound (ignoring delta)
############################################################

margin.bound <- function (m, d=1, error=1) {
  return (2/m * (64*d * log(exp(1)*m/(8*d),2) * log(32*m,2) + log(m,2)) -
          error);
}

############################################################
## computes the derivative of the margin bound (for fixed d)
############################################################

margin.bound.derivative <- function (m, d=1) {
  return (2 * (960*d*log(2)^2 - 128*d*log(m)*log(2) - 192*d*log(2) +
               64*d*log(m) - 64*d*log(m)^2 + 320*d*log(d)*log (2) +
               64*d*log(d)*log(m) - log(m)*log(2) + 64*d -64*d*log(d) +
               log(2)) / (m^2 * log(2)^2));
}

############################################################
## computes the minimal value of m for a given d
## in the margin bound case
############################################################

margin.min.training.set <- function (d, error=1, verbosity=0) {
  min.m <- integer (0);
  for (i in 1:length (d)) {
    m <- 30000 * d [i];
    if (verbosity > 0) {
      cat ("Determining the minimal training set size for d=", d [i], "\n");
    }
    while (abs (margin.bound (m, d [i], error)) > 1e-4) {
      if (verbosity > 0) {
        cat ("m=", m, "\tbound=", margin.bound (m, d [i], error), "\n");
      }
      m <- m - margin.bound (m, d [i], error) /
        margin.bound.derivative (m, d [i]);
    }
    if (verbosity > 0) {
      cat ("convergence suceeded\n\n");
    }
    min.m <- c (min.m, floor (m));
  }
  return (min.m);
}

############################################################
## plots the minimal required training sample size for
## the margin bound
############################################################

margin.make.plot <- function (d, error=1) {
  m <- margin.min.training.set (d, error=error);

  plot (d, m, type="n",
        xlab="margin complexity",
        ylab="minimal training set size",
        cex.axis = 2.5, cex.lab=3.5);
  lines (d, m, lwd=4);
  points (d, m, pch=19, lwd=2, cex=1.5);
}

############################################################
## plots the minimal required training set size
############################################################

book <- function (output='SCREEN') {

  ############################################################
  ## plots the VC bound over the whole range
  if (output == 'PS') {
    postscript ("VC_bound.ps")
  }
  par (mai=c (1.25,1.25,0.25,0.25));
  x <- seq (0, 1, length=1000);
  plot (x, x*(1+log(2)-log(x)), type="n", 
        cex.lab=3.5, cex.axis=2.5,
        xlab=expression(frac(nu,m)),
        ylab="complexity term");
  lines (x, x*(1+log(2)-log(x)), lwd=4);
  lines (x, x, lwd=4, lty=2);
  legend (0.0, 1.7, 
          c (expression (paste (frac(nu,m),
              bgroup("(",paste (ln, bgroup ("(",frac(2*m,nu)+1,")")),")"))),
             expression (frac(nu,m))),
          text.width=0.3,
          cex=2.5, lty=c(1,2), lwd=4);
  readline ("Press any key for next plot");
  dev.off ();
  if (output == 'PS') {
    cat ("[VC_bound.ps generated]\n");
  }

  ############################################################
  ## plots the VC bound in the small complexity range
  if (output == 'PS') {
    postscript ("VC_bound2.ps");
  }
  par (mai=c (1.25,1.25,0.25,0.25));
  x <- seq (0, 1/30, length=1000);
  plot (x, x*(1+log(2)-log(x)), type="n", 
        cex.lab=3.5, cex.axis=2.5,
        xlab=expression(frac(nu,m)),
        ylab="complexity term");
  lines (x, x*(1+log(2)-log(x)), lwd=4);
  lines (x, 5*x, lwd=4, lty=2);
  legend (0.0, 0.165, 
          c (expression (paste (frac(nu,m),
              bgroup("(",paste (ln, bgroup ("(",frac(2*m,nu)+1,")")),")"))),
             expression (5 * frac(nu,m))),
          text.width=0.01,
          cex=2.5, lty=c(1,2), lwd=4);
  readline ("Press any key for next plot");
  dev.off ();
  if (output == 'PS') {
    cat ("[VC_bound2.ps generated]\n");
  }

  ############################################################
  ## plots the fat shattering bound of minimal training sample size
  if (output == 'PS') {
    postscript (file="fat_bound_minimal_sample_size.ps");
  }
  par (mai=c (1.15,1.25,0.25,0.25));
  d <- seq (1, 10, by=1);
  fat.make.plot (d)
  readline ("Press any key for next plot");
  dev.off ();
  if (output == 'PS') {
    cat ("[fat_bound_minimal_sample_size.ps generated]\n");
  }

  ############################################################
  ## plots the margin shattering bound of minimal training sample size
  if (output == 'PS') {
    postscript (file="margin_bound_minimal_sample_size.ps");
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  d <- seq (1, 10, by=1);
  margin.make.plot (d)
  readline ("Press any key for next plot");
  dev.off ();
  if (output == 'PS') {
    cat ("[margin_bound_minimal_sample_size.ps generated]\n");  
  }
}

book()