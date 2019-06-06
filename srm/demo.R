###
### demo        creates all plots for the structural risk minimisation
###             section
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## BOOK plots
############################################################

book <- function (output='SCREEN') {
  if (output == 'PS') {
    postscript (file="srm.ps");
  }
  par (mai=c (0.75,1.0,0.25,0.25));

  ## global parameters  
  delta <- 0.05;
  m <- 5000;
  d <- seq (1, 20, by=1);

  ## compute VC term
  VC <- sqrt (8/m * (log (4 / (delta / 20)) + d * log (2*exp(1)*m/d)));

  ## training error sequence
  error <- round (1.5^(-d)* 1.5 * 1/2 * 5000);

  ## compute sum
  bound <- error/m + VC;

  ## plot the results
  plot (c (0, 20),
        c (0, max (bound)),
        type="n", cex.lab=2.5, cex.axis=2.0,
        xlab="model index", ylab="training/generalization error (bound)")

  ## draw the three lines
  lines (d, error/m, lwd=2, lty=1);
  lines (d, VC, lwd=2, lty=1);
  lines (d, bound, lwd=6, lty=1);
  
  points (d, error/m, pch=19, cex=2);
  points (d, VC, pch=19, cex=2);
  points (d, bound, pch=19, cex=2);

  ## mark the best model
  idx <- order (bound);
  points (idx [1], bound [idx [1]], pch=4, cex=3, lwd=4);

  x <- seq (-1, idx [1], length=100);
  lines (x, rep (bound [idx [1]], length=100), lty=3); 

  y <- seq (-1, bound [idx [1]], length=100);
  lines (rep (idx [1], length=100), y, lty=3); 

  ## label the curves
  text (1.5, bound [1], "bound", cex=1.8, adj=0);
  text (9, VC [8] - 0.01, "VC complexity term", cex=1.8, adj=0);
  text (11, 0.02 + error [10]/m, "training error", cex=1.8, adj=0);
  
  readline ("Press any key for next plot");
  dev.off ();
}

book ();
