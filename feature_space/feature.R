###
### feature.R           depicts mapping in feature space
###
### 2000 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

book <- function () {

  ## this code creates the plot of the embedding of
  ## the unite square in R^2 
  N <- 20;
  
  x1 <- seq (0, 1, length=N);
  x2 <- seq (0, 1, length=N);

  X1 <- x1;
  X2 <- x2^2;
  X3 <- matrix (0, nrow=N, ncol=N);

  for (i in 1:N) {
    for (j in 1:N) {
      X3 [i,j] <- x1 [i] * x2 [j];
    }
  }

  postscript (file="../../ps/feature.ps", horizontal=FALSE);
  par (mai=c (0,0.1,0,0));
  palette (gray (0:256/256));
  persp (X1, X2, X3,
         box=TRUE, border=gray(0),
         scale=TRUE, shade=0.7, 
         theta=30, phi=30,
         zlim=c(-0.5, 1.5),
         xlab="feature 1", ylab="feature 2", 
         zlab="feature 3", ticktype="detailed", cex=1.2);
  readline ("Press any key for next plot");
  dev.off ();
  
  ## this code creates possible decision surfaces in the original input
  ## space

  N  <- 1000;
  X1 <- seq (-1, 1, length=N);
  W1 <- c (-1, -1, -1,  0,  0,  0, +1, +1, +1);
  W3 <- c (-1,  0, +1, -1,  0, +1, -1,  0, +1);
  n  <- length (W3);
  W2 <- rep (1, length=n);
  X2.plus  <- matrix (0, nrow=n, ncol=N);
  X2.minus <- matrix (0, nrow=n, ncol=N);

  for (i in 1:n) {
    X2.plus [i,]  <- -W3 [i]/(2 * W2 [i]) * X1 +
      sqrt (X1 * (W3 [i]^2 * X1 - 4 * W1 [i] * W2 [i])/ (2 * W2 [i]));
    X2.minus [i,] <- -W3 [i]/(2 * W2 [i]) * X1 -
      sqrt (X1 * (W3 [i]^2 * X1 - 4 * W1 [i] * W2 [i])/ (2 * W2 [i]));
  }

  xr <- range (X1);
  x2 <- c (as.vector (X2.plus), as.vector (X2.minus));
  yr <- range (x2 [!is.nan (x2)]);
 
  postscript (file="../../ps/decision_surfaces.ps", horizontal=FALSE);
  par (mai=c (0.85,1.0,0.1,0.1));
  plot (xr, yr, type="n", xlab=expression (x[1]), ylab=expression (x[2]),
        cex.lab=3, cex.axis=2);

  types <- c (1, 2, 4);
  for (i in 1:n) {
    lines (X1, X2.plus  [i,], lwd=2, lty=types [W1[i]+2]);
    lines (X1, X2.minus [i,], lwd=2, lty=types [W1[i]+2]);
  }
  
  readline ("Press any key for next plot");
  dev.off ();
}

book ();


