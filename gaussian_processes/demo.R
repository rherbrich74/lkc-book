### demo        demonstrates Gaussian process regression
###             and evidence maximisation on some toy data
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
## load the plotting routines from rbf_demo
############################################################

source ("rbf_demo.R")

############################################################
## this routine computes the prediction + variance of a
## Gaussian process at a given set "test" of points
############################################################

GP.predict <- function (X,          ## data matrix of training points
                        y,          ## target values at the training points
                        test,       ## data matrix of test points
                        var=1,      ## assumed variance level
                        type="rbf", ## kernel type
                        ...) {      ## further kernel parameters

  ## parameter check
  if (is.null (dim (X)))
    stop ("Data matrix is not a matrix.");

  if (is.null (dim (test)))
    stop ("Data matrix is not a matrix.");

  if (nrow (X) != length (y))
    stop ("Number of training points and target values do not match.");

  if (ncol (X) != ncol (test))
    stop ("Number of features of training and test data do not match");

  if (var < 0)
    stop ("Variance has to be non-negative");

  ## compute the Gram matrix
  G <- kernel.linear (X, t(X), type=type, ...);

  ## compute the matrix of inner product between test points and
  ## training points
  H <- kernel.linear (test, t(X), type=type, ...);
  H2 <- kernel.linear (test, t(test), type=type, ...);
  
  ## invert the (Gram matrix + var*identiy);
  G2 <- G + diag (rep (var, length=nrow (X)));
  inv.G2 <- solve (G2);

  ## use the inverse to compute the predicitions + variances
  y.pred <- H %*% inv.G2 %*% cbind (y);
  v <- diag (H2) + var - diag (H %*% inv.G2 %*% t(H));  

  ## return the GP predictions
  return (list (y=y.pred, var=v));
}

############################################################
## this routine computes the log-evidence for a given training
## set and kernel + variance 
############################################################

GP.evidence <- function (X,          ## data matrix of training points
                         y,          ## target values at the training points
                         var=1,      ## assumed variance level
                         type="rbf", ## kernel type
                         G=NULL,     ## the full Gram matrix can be supplied
                                     ## if already computed
                         ...) {      ## further kernel parameters

  ## parameter check
  if (is.null (dim (X)))
    stop ("Data matrix is not a matrix.");

  if (nrow (X) != length (y))
    stop ("Number of training points and target values do not match.");

  if (var < 0)
    stop ("Variance has to be non-negative");

  ## compute the (Gram matrix + var*identiy), if necessary
  if (is.null (G)) {
    G <- kernel.linear (X, t(X), type=type, ...);
    G <- G + diag (rep (var, length=nrow (X)));
  }

  ## use the eigenvalue decompisition to speed up computations
  E <- eigen (G, symmetric=TRUE);
  inv.G <- E$vectors %*% diag (1/E$values) %*% t(E$vectors);
  det.G <- prod (E$values);

  evidence <- 1/sqrt (det.G) * exp (-1/2 * rbind (y) %*% inv.G %*% cbind (y));
  
  ## return the evidence
  return (evidence);
}

############################################################
## plots the prediction of a Gaussian process together with 
## errorbars on the outputs
############################################################

plot.GP <- function (X, y, type="rbf", var=0,
                     plot.data=TRUE, plot.var=TRUE, ...) {
  
  ## parameter check
  if (is.null (dim (X)))
    stop ("Data matrix is not a matrix.");

  if (nrow (X) != length (y))
    stop ("Number of training points and target values do not match.");

  if (ncol (X) > 2)
    stop ("Plots only work for one- and two-dimensional data.");

  if (var < 0)
    stop ("Variance has to be non-negative");

  #################### 1-D case ####################
  ## compute the GP-predictions and plot them
  if (ncol (X) == 1) {
    x.min <- min (X [,1]);
    x.max <- max (X [,1]);
    x <- cbind (seq (x.min - (x.max - x.min) * 0.3,
                     x.max + (x.max - x.min) * 0.3, length=100));
    res <- GP.predict (X, y, x, var=var, type=type, ...);
    
    plot (x, res$y,
          ylim=range (c (res$y - sqrt (res$var), res$y + sqrt (res$var))),
          type="l", lwd=3,
          xlab="x", ylab="t(x)", cex.lab=4.0, cex.axis=2.0);
    if (plot.data) {
      points (X, y, pch=19, cex=2);
    }
    if (plot.var) {
      lines (x, res$y - sqrt (res$var), lty=2);
      lines (x, res$y + sqrt (res$var), lty=2);
    }
  }

  #################### 2-D case ####################
  ## compute the GP-predictions and plot them
  if (ncol (X) == 2) {
    stop ("Not yet implemented.");
  }
  
  return;
}

############################################################
## EXAMPLE code - this code allows you to
##    1) input data point using the mouse pointer
##    2) computes the evidence
##    3) plots the local maxima solutions 
############################################################

evidence <- function (n=10,                 ## number of training points
                      N=20,                 ## discretisations grid
                      var=c(0,1),           ## range of variance
                      sigma=c(0.1,1),       ## range of bandwidth
                      var.plot=c(0,1),      ## list of variance values for
                      sigma.plot=c(0.1,1),  ## list of bandwidths values for
                      pch=c(19,4),          ## the charatcers to be used
                      file=NULL) {          
  
  ## gather the data points
  if (!file.exists ("train.dat")) {
    plot (c(0,6), c(-1, 3), type="n");
    axis (1);
    axis (2);
    D <- locator (n=n, type="p", pch=19);
    X <- cbind (D$x);
    y <- D$y;
    save (X, y, file="train.dat");
  } else {
    load (file="train.dat");
  }
  
  ## compute the evidence 
  var.list <- seq (var [1], var [2], length=N);
  sigma.list <- seq (sigma [1], sigma [2], length=N);
  evidence <- matrix (0, nrow=N, ncol=N);
  for (i in 1:N) {
    for (j in 1:N) {
      evidence [i,j] <- GP.evidence (X, y,
                                     var=var.list [i],
                                     sigma=sigma.list [j]);
    }
  }
  
  ## plot the evidence
  if (!is.null (file)) {
    postscript (file=file [1]);
    par (mai=c (1.0,1.25,0.25,0.25));
  } else {
    x11 ();
  }
  image (var.list, sigma.list, log (evidence+1e-4),
         col=gray (exp (seq (0,1,length=200))/exp(1)),
         xlab="variance", ylab="bandwidth",
         cex.axis=2, cex.lab=4);
  contour (var.list, sigma.list, log (evidence+1e-4),
           labcex=2, lwd=1.2,
           levels=c(-5, -4, -3, -2.5, -2, -1.5, -1.2),
           add=TRUE);

  ## mark the models specified
  for (i in 1:length (var.plot)) {
    points (var.plot [i], sigma.plot [i], pch=pch [i], cex=4, lwd=4);
  }
  readline ("Evidence plotted");
  dev.off ();

  ## cycle through all data plots
  for (i in 1:length (var.plot)) {
    if (!is.null (file)) {
      postscript (file=file [1+i]);
      par (mai=c (1.0,1.25,0.25,0.25));
    } else {
      x11 ();
    }
    plot.GP (X, y, var=var.plot [i], sigma=sigma.plot [i]);
    readline ("Data fit plot");
    dev.off ();
  }
}

############################################################
## returns a sample from a multidimensional Gaussians
## (taken from S-PLUS help system)
############################################################

rmultnorm <- function (n, mu, vmat) {
  p <- ncol(vmat);
  ans <- matrix (rnorm (n * p), nrow = n) %*% vmat;
  ans <- sweep (ans, 2, mu, "+");
  return (drop (ans));
}                                                                                 
############################################################
## samples from a Gaussian process 
############################################################

GP.sample <- function (N=20, sigma1=1, sigma2=1) {

  ## generate a grid
  x <- seq (-10, 10, length=N);
  y <- seq (-10, 10, length=N);

  ## generate the data matrix for the GP process
  X <- rep (x, rep (N, length=N));
  Y <- rep (y, length=N*N);
  data <- cbind (X, Y);

  G <- data %*% rbind (c(1/sigma1, 0), c(0, 1/sigma2)) %*% t(data);
  
  ## compute the covariance matrix
  C <- kernel (G, diag (G), diag (G), type="rbf", sigma=1/sqrt (2));

  ## sample the Gaussian
  z <- matrix (rmultnorm (1, 0, C), nrow=N, ncol=N);
  
  ## return the sample
  return (list (x=x, y=y, z=z));
}

############################################################
## BOOK plots
############################################################

book <- function () {
  ############################################################
  ## generate the index plots 
  if (!file.exists ("../../ps/GP_evidence.ps")) {
    evidence (var=c (0, 0.6), sigma=c (0.1, 3.5), var.plot=c (0, 0.5),
              sigma.plot=c(1.1, 3),
              file=c("../../ps/GP_evidence.ps",
                "../../ps/GP_nonoise_model.ps",
                "../../ps/GP_fullnoise_model.ps"), N=50);
  }

  ############################################################
  ## plot samples from the ARD prior 
  if (!file.exists ("../../ps/ARD_sample1.ps")) {
    set.seed (0);
    postscript (file="../../ps/ARD_sample1.ps");
    persp (GP.sample (N=40, sigma1=5, sigma2=5),
           theta=35, phi=30, shade=0.7, border=gray (0.25),
           ticktype="detailed", xlab="input dimension 1", ylab="input dimension 2",
           zlab="function values", cex=1.3,
           lphi=30, ltheta=-85);
    cat ("[ARD_sample1.ps plotted]\n");
    dev.off ();
    postscript (file="../../ps/ARD_sample2.ps");
    persp (GP.sample (N=40, sigma1=5, sigma2=100),
           theta=35, phi=30, shade=0.7, border=gray (0.25),
           ticktype="detailed", xlab="input dimension 1", ylab="input dimension 2",
           zlab="function values", cex=1.3,
           lphi=30, ltheta=-85);
    cat ("[ARD_sample2.ps plotted]\n");
    dev.off ();
  }
  
  ############################################################
  ## plot the effect of a latent variable
  if (!file.exists ("../../ps/input_signal.ps")) {

    ## plot the input signal
    postscript (file="../../ps/input_signal.ps", horizontal=FALSE);
    par (mai=c (1.0,1.25,0.25,0.25));
    load ("train.dat");
    x.min <- min (X [,1]);
    x.max <- max (X [,1]);
    x <- cbind (seq (x.min - (x.max - x.min) * 0.3,
                     x.max + (x.max - x.min) * 0.3, length=100));
    out <- GP.predict (X, y, x, var=0)$y;
    plot (x, out, type="l", lwd=3, xlab="x", ylab="t(x)",
          cex.lab=4.0, cex.axis=2.5);
    readline ("[../../ps/input_signal.ps plotted]\n");
    dev.off ();

    ## plot different sigmoid's
    postscript (file="../../ps/sigmoid.ps", horizontal=FALSE);
    par (mai=c (1.0,1.25,0.25,0.25));
    beta.list <- c (0.1, 1, 5);
    z <- seq (-max (abs (y)), max (abs (y)), length=100);
    plot (range (z), c (0, 1), type="n", xlab="t", ylab=expression (pi(t)),
          cex.lab=4, cex.axis=2.5);
    for (i in 1:length (beta.list)) {
      lines (z, exp (2/beta.list [i] * z) / (1 + exp (2/beta.list [i] * z)),
             lty=i, lwd=3);
    }
    legend (x [1], 1.0, c (expression (paste (beta, "=0.1")),
                           expression (paste (beta, "=1.0")),
                           expression (paste (beta, "=5.0"))),
            lty=1:3, lwd=3, cex=3);
    readline ("[../../ps/sigmoid.ps plotted]\n");
    dev.off ();

    ## plot the transfered latent functions
    postscript (file="../../ps/transfered_input_signal.ps", horizontal=FALSE);
    par (mai=c (1.0,1.25,0.25,0.25));
    plot (range (x), c (0, 1), type="n", xlab="x",
          ylab=expression (paste (bold(P)[paste (Y, group ("|", paste ("X=x"), ""))],"(+1)")),
          cex.lab=4, cex.axis=2.5);
    for (i in 1:length (beta.list)) {
      out <- GP.predict (X, y, cbind (x), sigma=0.5, var=0)$y;
      lines (x,
             exp (2/beta.list [i] * out) / (1 + exp (2/beta.list [i] * out)),
             lty=i, lwd=3);
    }
    readline ("[../../ps/transfered_input_signal.ps plotted]\n");
    dev.off ();
  }
}

