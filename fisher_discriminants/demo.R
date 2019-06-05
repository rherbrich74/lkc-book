###
### demo        creates all plots for the fisher discriminant section
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## main parameters
############################################################

D <- 4.5;                      # dimension of output window (generative plot)

D2.x <- 4.3;                   # dimension of output window (projections plot)
D2.y <- 5.8;                   

N <- 40;                       # number of samples per classa

N2 <- 10;                      # number of samples per classa

mu <- cbind (c (-1, +1),       # mean vectors
             c (+1, -1));
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
  
  ############################################################
  ## generate "projection" prior plots
  if (!file.exists ("../../ps/FD_generative.ps")) {
    postscript ("../../ps/FD_generative.ps", horizontal=FALSE);
    
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
    
    ## compute the Fisher discriminant
    w <- Sigma.inv %*% (mu.plus - mu.minus);
    b <- as.vector (1/2 * (t (mu.minus) %*% Sigma.inv %*% mu.minus -
                           t (mu.plus) %*% Sigma.inv %*% mu.plus));
    
    ## compute the denominator of the Gaussian density
    g.norm <- 1/((2 * pi) * sqrt (det (Sigma))); 
    
    ## create the plot window
    plot.new ();
    plot.window (c (-D, D), c(-D, D));

    ## draw the dataset on the screen
    points (X.plus [, 1], X.plus [, 2], pch=20, cex=2, col=gray (0.5));
    points (X.minus [, 1], X.minus [, 2], pch=4, lwd=2, cex=2, col=gray (0.5));

    ## draw the decision surface
    x <- seq (-D, D, length=1000);
    y <- -(w [1] * x + b) / w [2];
    lines (x, y, lwd=4);

    ## draw the means
    points (mu.plus [1], mu.plus [2],   pch=17, cex=2, lwd=2);
    points (mu.minus [1], mu.minus [2], pch=17, cex=2, lwd=2);

    ## draw the concentration ellipses
    x <- seq (-D, D, length=100);
    y <- seq (-D, D, length=100);
    z.plus <- outer (x, y, function (x, y) {
      D <- (x - mu.plus [1])^2 * Sigma.inv [1, 1] +
        (y - mu.plus [2])^2 * Sigma.inv [2, 2] +
          2 * (x - mu.plus [1]) * (y - mu.plus [2]) * Sigma.inv [1, 2];
      return (g.norm * exp (-1/2 * D));
    });
    z.minus <- outer (x, y, function (x, y) {
      D <- (x - mu.minus [1])^2 * Sigma.inv [1, 1] +
        (y - mu.minus [2])^2 * Sigma.inv [2, 2] +
          2 * (x - mu.minus [1]) * (y - mu.minus [2]) * Sigma.inv [1, 2];
      return (g.norm * exp (-1/2 * D));
    });
    
    contour (x, y, z.plus,  levels=c (0.0028), drawlabels=FALSE,
             lwd=2, add=TRUE);
    contour (x, y, z.minus, levels=c (0.0028), drawlabels=FALSE,
             lwd=2, add=TRUE);

    ## draw the ellipses that touch each other
    contour (x, y, z.plus,  levels=c (0.00272), drawlabels=FALSE,
             lwd=1, lty=2, add=TRUE);
    contour (x, y, z.minus, levels=c (0.00272), drawlabels=FALSE,
             lwd=1, lty=2, add=TRUE);
    
    dev.off ();
  }

  ############################################################
  ## generate "data distribution" prior plots
  if (!file.exists ("../../ps/FD_projections.ps")) {
    postscript ("../../ps/FD_projections.ps", horizontal=FALSE);
    
    ## create the datasets
    X.plus <-  rmultnorm (N2, mu [, 1], S.plus);
    X.minus <- rmultnorm (N2, mu [, 2], S.minus);
    X <- rbind (X.plus, X.minus);
    
    ## compute the means
    mu.plus  <- 1/N2 * apply (X.plus,  2, sum);
    mu.minus <- 1/N2 * apply (X.minus, 2, sum);
    
    ## compute the covariance
    Sigma <- t (X) %*% X - N2 * mu.plus %*% t (mu.plus) -
      N2 * mu.minus %*% t (mu.minus);
    Sigma.inv <- solve (Sigma);
    
    ## compute the Fisher discriminant
    w <- Sigma.inv %*% (mu.plus - mu.minus);
    b <- as.vector (1/2 * (t (mu.minus) %*% Sigma.inv %*% mu.minus -
                           t (mu.plus) %*% Sigma.inv %*% mu.plus));
    
    ## compute the denominator of the Gaussian density
    g.norm <- 1/((2 * pi) * sqrt (det (Sigma))); 
    
    ## create the plot window
    plot.new ();
    plot.window (c (-D2.x, D2.x), c(-D2.y, D2.y));

    ## draw the dataset on the screen
    points (X.plus [, 1], X.plus [, 2], pch=20, cex=2, col=gray (0.4));
    points (X.minus [, 1], X.minus [, 2], pch=3, lwd=2, cex=2, col=gray (0.4));

    ## draw the decision surface
    x <- seq (-D2.x, D2.x, length=1000);
    y <- -(w [1] * x + b) / w [2];
    lines (x, y, lwd=4);

    ## compute the offset vector
    ofs      <- 4 *  cbind (-w [2], w [1]);
    ofs.mean <- 5 * cbind (-w [2], w [1]);
    ofs.var  <- 6 * cbind (-w [2], w [1]);
    
    ## draw the orthogonal line
    y <- w [2] / w [1] * x;
    x <- x + ofs [1]; y <- y + ofs [2];
    lines (x, y, lwd=2);

    ## project the points
    w.norm <- as.vector (t (w) %*% w);
    l.plus <- (X.plus %*% w) / w.norm;
    l.minus <- (X.minus %*% w) / w.norm;

    X.plus.proj  <- cbind (w [1] * l.plus + ofs [1],
                           w [2] * l.plus + ofs [2]);
    X.minus.proj <- cbind (w [1] * l.minus + ofs [1],
                           w [2] * l.minus + ofs [2]);
    
    ## draw the projected points
    points (X.plus.proj  [, 1], X.plus.proj  [, 2] ,  pch=20, cex=2);
    points (X.minus.proj [, 1], X.minus.proj [, 2] ,  pch=3, lwd=2, cex=2);

    ## draw the projection lines
    for (i in 1:N2) {
      arrows (X.plus [i, 1], X.plus [i, 2],
              X.plus.proj [i, 1], X.plus.proj [i, 2],
              angle=10, lwd=1.5, length=0.1, col=gray (0.5));

      arrows (X.minus [i, 1], X.minus [i, 2],
              X.minus.proj [i, 1], X.minus.proj [i, 2],
              angle=10, lwd=1.5, length=0.1, col=gray (0.5));
    }

    ## compute the mean and variance of the projection
    mu.plus  <- mean (l.plus);
    mu.minus <- mean (l.minus);

    mu.plus.proj <- rbind (w [1] * mu.plus + ofs.mean [1],
                           w [2] * mu.plus + ofs.mean [2]);
    
    mu.minus.proj <- rbind (w [1] * mu.minus + ofs.mean [1],
                            w [2] * mu.minus + ofs.mean [2]);

    ## draw the mean projection
    points (mu.plus.proj  [1], mu.plus.proj  [2], pch=17, cex=2, lwd=2);
    points (mu.minus.proj [1], mu.minus.proj [2], pch=17, cex=2, lwd=2);
    arrows (mu.plus.proj  [1], mu.plus.proj  [2],
            mu.minus.proj [1], mu.minus.proj [2],
            code=3, lwd=2, angle=15, length=0.1, col=gray (0.5));
    
    ## compute the mean and variance of the projection
    std.plus  <- sqrt (var (l.plus));
    std.minus <- sqrt (var (l.minus));

    std.plus.left  <- rbind (w [1] * (mu.plus - std.plus) + ofs.var [1],
                             w [2] * (mu.plus - std.plus) + ofs.var [2]);    
    std.plus.right <- rbind (w [1] * (mu.plus + std.plus) + ofs.var [1],
                             w [2] * (mu.plus + std.plus) + ofs.var [2]);
    
    std.minus.left  <- rbind (w [1] * (mu.minus - std.minus) + ofs.var [1],
                              w [2] * (mu.minus - std.minus) + ofs.var [2]);
    std.minus.right <- rbind (w [1] * (mu.minus + std.minus) + ofs.var [1],
                              w [2] * (mu.minus + std.minus) + ofs.var [2]);
    
    ## draw the mean projection
    arrows (std.plus.left  [1], std.plus.left  [2],
            std.plus.right [1], std.plus.right [2],
            code=3, lwd=2, angle=15, length=0.1, col=gray (0.5));

    arrows (std.minus.left  [1], std.minus.left  [2],
            std.minus.right [1], std.minus.right [2],
            code=3, lwd=2, angle=15, length=0.1, col=gray (0.5));
    
    dev.off ();
  }

}

book ()
