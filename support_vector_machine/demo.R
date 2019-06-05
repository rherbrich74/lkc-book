###
### demo        creates all plots for the support vector machine
###             section (must be used in conjunction with xfig)
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## main parameters
############################################################

D.x <- 6.5;                    # dimension of output window (generative plot)
D.y <- 8.9;                    

alpha <- 35 * pi/180;          # rotation angle of the second hyperplane (rad)

N <- 40;                       # number of samples per classa

mu <- cbind (c (-3, +3),       # mean vectors
             c (+3, -3));
C <- list (rho=c (0.5, 0.5),   # covariance matrices
          sigma.x=c (1, 1),
           sigma.y=c (1, 1));


############################################################
## a function to generate random multivariate Gaussians
## directly taken from S-PLUS help files
############################################################
rmultnorm <- function (n, mu, vmat) {
  p <- ncol (vmat);
  ans <- matrix (rnorm (n * p), nrow = n) %*% vmat;
  ans <- sweep (ans, 2, mu, "+");
  drop (ans);
}


############################################################
## this function draw the margin plot
############################################################
draw.plane <- function (X.plus, X.minus, w) {
  ## create the plot window
  plot.new ();
  plot.window (c (-D.x, D.x), c(-D.y, D.y));
  
  ## determine the closest data points
  out <- rbind (X.plus %*% w, -X.minus %*% w);
  gamma <- min (out);
  
  ## compute the corners to fill the dead zone
  x <- double (0);
  y <- double (0);
  
  x <- c (y, -(2 * D.y * w [2] - gamma)/ w [1]);
  y <- c (y, 2 * D.y);
  
  x <- c (x, 2 * D.x);
  y <- c (y, 2 * D.y);
  
  x <- c (x, 2 * D.x);
  y <- c (y, -(2 * D.x * w [1] + gamma)/ w [2]);
  
  x <- c (x, -(-2 * D.y * w [2] + gamma)/ w [1]);
  y <- c (y, -2 * D.y);
  
  x <- c (x, -2 * D.x);
  y <- c (y, -2 * D.y);
  
  x <- c (x, -2 * D.x);
  y <- c (y, -(-2 * D.x * w [1] - gamma)/ w [2]);
  
  ## draw the dead zone
  polygon (x, y, border=0, col=gray (0.9));
  
  ## draw the dataset on the screen
  points (X.plus [, 1], X.plus [, 2], pch=4, lwd=2, cex=2, col=gray (0.5));
  points (X.minus [, 1], X.minus [, 2], pch=20, cex=2, col=gray (0.5));
  
  ## draw the decision surface
  x <- seq (-2 * D.x, 2 * D.y, length=1000);
  y <- -(w [1] * x) / w [2];
  lines (x, y, lwd=4);
  
  ## draw the margin lines
  y <- -(w [1] * x + gamma) / w [2];
  lines (x, y, lty=2, lwd=2);
  
  y <- -(w [1] * x - gamma) / w [2];
  lines (x, y, lty=2, lwd=2);

  ## compute the offset vectors
  ofs <- c (-w [2], w [1]);
  
  ## draw the weight vector
  arrows (0, 0, w [1] * gamma, w [2] * gamma, lwd=2,
          angle=10, length=0.15, code=3);
  arrows (ofs [1], ofs [2],
          ofs [1] + 2 * w [1] / gamma, ofs [2] + 2 * w [2] / gamma,
          lwd=2, angle=10, length=0.15);
  
}
############################################################
## BOOK plots
############################################################

book <- function () {

  ## fix the random number generator to reproduce "random" numbers
  set.seed (4);
  
  ## generate the covariance matrices
  S.plus <- diag (c (C$sigma.x [1], C$sigma.y [1])) %*%
    cbind (c (1, C$rho [1]), c (C$rho [1], 1)) %*%
      diag (c (C$sigma.x [1], C$sigma.y [1]));
  
  S.minus <- diag (c (C$sigma.x [2], C$sigma.y [2])) %*%
    cbind (c (1, C$rho [2]), c (C$rho [2], 1)) %*%
      diag (c (C$sigma.x [2], C$sigma.y [2]));
  
  ## create the datasets
  X.plus <-  rmultnorm (N, mu [, 1], S.plus);
  X.minus <- rmultnorm (N, mu [, 2], S.minus);
  X <- rbind (X.plus, X.minus);
  
  ## compute the means
  mu.plus  <- 1/N * apply (X.plus,  2, sum);
  mu.minus <- 1/N * apply (X.minus, 2, sum);
  
  ## compute the covariance
  Sigma <- t (X) %*% X - N * mu.plus %*% t (mu.plus) -
    N * mu.minus %*% t (mu.minus);
  Sigma.inv <- solve (Sigma);
  
  ############################################################
  ## generate hyperplane plots
  if (!file.exists ("../../ps/svm_hyperplane1.ps")) {
    postscript ("../../ps/svm_hyperplane1.ps", horizontal=FALSE);
     
    ## compute and draw discriminant no 1
    w <- Sigma.inv %*% (mu.plus - mu.minus);
    w <- 1 / as.vector (sqrt (t (w) %*% w)) * w;

    draw.plane (X.plus, X.minus, w);
    
    dev.off ();
  }

  if (!file.exists ("../../ps/svm_hyperplane2.ps")) {
    postscript ("../../ps/svm_hyperplane2.ps", horizontal=FALSE);
     
    ## compute and draw discriminant no 2
    w <- cbind (c (cos (alpha), -sin (alpha)),
                c (sin (alpha), cos (alpha))) %*% w;

    draw.plane (X.plus, X.minus, w);
    
    dev.off ();
  }

}






