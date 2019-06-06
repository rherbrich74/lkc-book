### demo        creates all plots for the relevance vector
###             machines section
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## computes the 'marginalised' weight vector prior
############################################################

marginalised.log.prior <- function (x, a=0, b=1e10) {

  ## parameter checks
  if (!is.null (dim (x))) {
    stop ("Only vector of number allowed as the first argument.");
  }

  if (a < 0 || b < 0) {
    stop ("Gamma prior parameters out of range.");
  }

  ## compute the function
  return (lgamma (a + 1/2) -  log (sqrt (2 * pi)) - lgamma (a) - a * log (b) 
          -(a + 1/2) * log (1 / b + x^2 / 2));
}

############################################################
## BOOK plots
############################################################

book <- function (output='SCREEN') {
  
  ############################################################
  ## generate "marginalised" prior plots
  if (output == 'PS') {
    postscript ("rvm_marginalised_prior.ps");
  }
  x <- seq (-2, 2, length=1000);
  
  m1 <- marginalised.log.prior (x, a=0.01, b=100);
  m2 <- marginalised.log.prior (x, a=0.001, b=1000);
  m3 <- marginalised.log.prior (x, a=0.0001, b=10000);
  

  par (mai=c (1.0,1.25,0.25,0.25));
  xlim <- range (x);
  ylim <- range (c (m1, m2, m3));
  plot (xlim, ylim, type="n", xlab="w",
        ylab=expression (paste (log, group ("(", paste (bold (f) [W[i]],
            group ("(", "w", ")")), ")"))),
        cex.lab=3.5, cex.axis=2.5);
  lines (x, m1, lwd=4, lty=3, col=gray (0.5));
  lines (x, m2, lwd=4, lty=2, col=gray (0.3));
  lines (x, m3, lwd=4, lty=1, col=gray (0.1));
  legend (0.5, ylim [2] - 0.1, c ("a=1e-2, b=1e2",
                                  "a=1e-3, b=1e3",
                                  "a=1e-4, b=1e4"),
          lwd=3, lty=c (3, 2, 1),
          col=c(gray (0.5), gray (0.3), gray (0.1)),
          cex=2.0);
  
  if (output == 'PS') {
    cat ("rvm_marginalised_prior.ps created\n");
  } else {
    readline ("Press any key to continue");
  }
  dev.off ();

  ############################################################
  ## generate "marginalised" 2D priors
  if (output == 'PS') {
    postscript ("rvm_marginalised_prior2.ps");
  }
  
  x <- seq (-2, 2, length=50);
  y <- seq (-2, 2, length=50);
  d <- outer (x, y, function (a, b) {
    return (marginalised.log.prior (a, a=0.00001, b=10000) +
            marginalised.log.prior (b, a=0.00001, b=10000)); });
  persp (x, y, d, r=1.2, phi=30, theta=135, shade=0.8,
         xlab="second component of w",
         ylab="first component of w",
         zlab="marginalised prior density",
       cex=1.4, ticktype="detailed");

  if (output == 'PS') {
    cat ("rvm_marginalised_prior2.ps created\n");
    cat ("rvm_marginalised_prior.ps created\n");
  } else {
    readline ("Press any key to continue");
  }
  dev.off ();

}

book()