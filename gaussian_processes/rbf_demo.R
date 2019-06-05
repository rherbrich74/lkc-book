### rbf_demo        demonstrates the ridge problem as well as
###                 the noise problem
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
## computes the inner product matrix between two sets of m
## and n points using kernels. Instead of assuming given matrices
## it is just assumed that we know the inner product of the
## m x n items as well as their squared norm. The kernel parameters
## are given as single arguments.
############################################################

kernel <- function (XY, X2, Y2, type="rbf", sigma=1, degree=1, c=0) {
  
  ## parameter check
  if (nrow (XY) != length (X2))
    stop ("Matrix dimensions of XY and X2 does not match.");

  if (ncol (XY) != length (Y2))
    stop ("Matrix dimension of XY and Y2 does not match.");

  if (type != "rbf" && type != "linear" && type != "polynomial") 
    stop ("Kernel types must match 'rbf', 'polynomial' or 'linear'.");
  
  ## compute the RBF kernel
  if (type == "rbf") {
    m <- nrow (XY);
    n <- ncol (XY);

    H1 <- matrix (X2, nrow=m, ncol=n);
    H2 <- matrix (Y2, nrow=m, ncol=n, byrow=TRUE);

    G <- exp (-(H1 -2 * XY + H2) / (2*sigma^2));
  }

  ## compute the polynomial kernel
  if (type == "polynomial") {
    G <- (XY + c)^degree;
  }

  ## compute the linear kernel
  if (type == "linear") {
    G <- XY;
  }
  
  ## return Gram matrix
  class (G) <- "gram";
  return (G);
}

############################################################
## computes the inner product matrix between two matrices using kernels.
## Each row of the matrix X and each column of the matrix Y are
## assumed to be data items in vectorial form.
############################################################

kernel.linear <- function (X, Y, ...) {

  ## parameter check
  if (ncol (X) != nrow (Y)) 
    stop ("Dimensions of matrices do not match.");

  return (kernel (X %*% Y, diag (X %*% t(X)),
                       diag (t(Y) %*% Y), ...));
}

############################################################
## normalises a given Gram matrix
############################################################

normalise <- function (G) {

  ## parameter check
  if (class (G) != "gram")
    stop ("Argument has to be of class 'gram'.");
  
  ## normalise kernel matrix
  n <- sqrt (1/diag (G));
  G2 <- G * (n %*% t (n));
  class (G2) <- "gram";
  
  return (G2);
}

############################################################
## generates a two class problem in 2D by two Gaussian at (-M,-M) and
## (+M,+M) with standard deviation sigma (or sigma.plus and sigma.minus
## for the two different classes). Each class contains m (or m.plus and
## m.minus) number of examples and are ordered consecutively.
############################################################

gen.data <- function (M = 5, sigma = 1, sigma.plus = sigma,
                      sigma.minus = sigma, m = 50, m.plus = m,
                      m.minus = m) {

  ## generate the x and y coordinates separately 
  x <- c (rnorm (m.plus,  mean=M,  sd=sigma.plus),
          rnorm (m.minus, mean=-M, sd=sigma.minus));
  y <- c (rnorm (m.plus,  mean=M,  sd=sigma.plus),
          rnorm (m.minus, mean=-M, sd=sigma.minus));

  ## concatenate them and return the matrix
  ret <- list (X = cbind (x, y),
               class = c (rep (+1, m.plus), rep (-1, m.minus)));
  class (ret) <- "toy.data"
  return (ret);
}

############################################################
## generates a two class problem in 2D by four Gaussian as
## (+/-M,+/-MM) with standard deviation sigma (or sigma.plus and
## sigma.minus for the two different classes). Each class contains m
## (or m.plus and m.minus) number of examples and are ordered
## consecutively.
############################################################

gen.parity.data <- function (M = 5, sigma = 1, sigma.plus = sigma,
                             sigma.minus = sigma, m = 50, m.plus = m,
                             m.minus = m) {

  ## generate the x and y coordinates separately 
  x <- c (rnorm (m.plus / 2,  mean=M,  sd=sigma.plus),
          rnorm (m.plus / 2,  mean=-M,  sd=sigma.plus),
          rnorm (m.minus / 2, mean=M, sd=sigma.minus),
          rnorm (m.plus / 2,  mean=-M, sd=sigma.plus));
  y <- c (rnorm (m.plus / 2,  mean=M,  sd=sigma.plus),
          rnorm (m.plus / 2,  mean=-M, sd=sigma.minus),
          rnorm (m.minus / 2, mean=-M,  sd=sigma.plus),
          rnorm (m.minus / 2, mean=M, sd=sigma.minus));

  ## concatenate them and return the matrix
  ret <- list (X = cbind (x, y),
               class = c (rep (+1, m.plus), rep (-1, m.minus)));
  class (ret) <- "toy.data"
  return (ret);
}

############################################################
## plots a circle with a prespecified radius and center
############################################################

read.dataset <- function (file) {
  
  ## parameter check
  if (!file.exists (file))
    stop ("File does not exists!");

  ## read the header
  header <- scan (file, n = 3, quiet=TRUE);
  if (header [3] != 0)
    stop ("File is of incorrect format!");
  d <- header [1];
  n <- header [2];

  ## read the body
  D <- matrix (scan (file, skip=1, n=(d+1)*n, quiet=TRUE),
               nrow=n, byrow=TRUE);
  return (list (X=D[,1:d], class=D[,(d+1)]));
}

############################################################
## plots a circle with a prespecified radius and center
############################################################

circle <- function (r=1, x=0, y=0, fill=FALSE, ...) {
  ## parameter check
  if (r < 0) 
    stop ("Radius must be positive.");

  ## generate the x,y coordinates
  i <- seq (0, 2*pi, length=100);
  X <- x + sin (i) * r;
  Y <- y + cos (i) * r;

  ## plot the circle
  if (fill == FALSE) {
    lines (X, Y, ...);
  } else {
    polygon (X, Y, ...);
  }
}

############################################################
## plots 2D data in an single window
############################################################

plot.toy.data <- function (data, sd=0, fill=TRUE, ...) {
  
  ## parameter checks
  if (class (data) != "toy.data") 
    stop ("Incorrect argument class.");

  ## check if an X11 terminal is already open?
  if (dev.cur () == 1) 
    plot.new ();

  ## setup the plot
  plot (range (data$X [,1] + sd, data$X [,1] - sd),
        range (data$X [,2] + sd, data$X [,2] - sd),
        type="n",
        xlab = expression (x[1]),
        ylab = expression (x[2]),
        cex.lab = 1.8, cex.axis = 1.6);

  ## determine the color (fill with different colors but draw with the
  ## same color
  col.plus <- gray (0.8);
  col.minus <- gray (0.4);
  if (fill == FALSE) {
    col.plus <- gray (0.4);
    col.minus <- gray (0.4);
  }
  
  ## plot the data (circles?)
  for (i in 1:dim (data$X) [1]) {
    if (data$class [i] == -1) {
      if (sd != 0) 
        circle (r=sd, x=data$X [i,1], y=data$X [i,2], fill=fill,
                col=col.minus, ...);
      points (data$X [i,1], data$X [i,2], pch=19, cex=1.2, lwd=2, ...);
    } else {
      if (sd != 0) 
        circle (r=sd, x=data$X [i,1], y=data$X [i,2], fill=fill,
                col=col.plus, ...);
      points (data$X [i,1], data$X [i,2], pch=4, cex=1.2, lwd=2, ...);
    }
  }
}

############################################################
## plots a gram matrix (data) either by a surface plot
## or as an intensity image; additionally it is possible to
## include a colorbar; if matlab=TRUE then an .m file is
## generated which directly plots the surface
############################################################

plot.gram <- function (G, surface=TRUE, colorbar=FALSE,
                       matlab=NA, file=NA, ...) {

  ## parameter check
  if (nrow (G) != ncol (G)) 
    stop ("Data has to be a symmetric matrix.");

  if (!is.na (matlab) && surface == FALSE)
    stop ("Only surface plots can be exported to MATLAB.");
    
  ## surface plot
  if (surface == TRUE) {
    m <- dim (G) [1];

    if (is.na (matlab)) {
      ## print the file, if necessary
      if (!is.na (file)) {
        postscript (file=file, horizontal=FALSE);
      } else { 
        par (fig=c (0.05, 1.0, 0.05, 1.0));
      }
      persp (x=1:m, y=1:m, z=G,
             theta=330, phi=25,
             ltheta=120, lphi=30,
             shade=0.75, border=NA,
             scale=TRUE, ticktype="detailed",
             xlab="row index", ylab="column index",
             zlab="similarity", cex=1.2, ...)
      
      ## print the file, if necessary
      if (!is.na (file)) {
        dev.off ();
      }
      
    } else {
      ## MATLAB export
      ##############################################################
      file.name <- matlab;

      ## export header
      cat ("%% This file is automatically constructed by plot.gram\n\n",
           file=file.name);

      ## export data
      cat ("X = [", file=file.name, append=TRUE);
      for (i in 1:m) {
        cat ("[", G [i,], "];\n", file=file.name, append=TRUE);
      }
      cat ("];\n", file=file.name, append=TRUE);

      ## plotting routines 
      cat ("surfl (X); colormap gray; shading interp; \n", file=file.name,
           append=TRUE);
      cat ("set (gca, 'FontSize', 12);\n", file=file.name, append=TRUE);
      cat ("h=xlabel ('row index'); set (h, 'FontSize', 14); clear h;\n",
           file=file.name, append=TRUE);
      cat ("h=ylabel ('column index'); set (h, 'FontSize', 14); clear h;\n",
           file=file.name, append=TRUE);
      cat ("h=zlabel ('similarity'); set (h, 'FontSize', 14); clear h;\n",
           file=file.name, append=TRUE);
      cat ("axis tight\n", file=file.name, append=TRUE);

      ## plot in file, if necessary
      if (!is.na (file)) {
        cat ("print -deps ", file, "\n", file=file.name, append=TRUE);
        cat ("close all\n", file=file.name, append=TRUE);
      }
    }
  }
  
  ## image plot
  if (surface == FALSE) {
    m <- dim (G) [1];
    
    ## print the file, if necessary
    if (!is.na (file)) {
      postscript (file=file, horizontal=FALSE);
    } else {
      par (fig=c (0.05, 1.0, 0.05, 1.0));
    }
    
    image (1:m, 1:m, G,
           xlab="row index", ylab="column index",
           cex.lab=1.8, cex.axis=1.6,
           col=gray ((0:256)/256), ...);
    box ();
    
    ## print the file, if necessary
    if (!is.na (file)) {
      dev.off ();
    }
  }
  
}

############################################################
## EXAMPLE code
############################################################

exmp <- function (M=5, m=50, sigma=3, type="rbf", surface=TRUE,
                  matlab=NA, file=NA) {
  D <- gen.parity.data (M=M, m=m);
  if (type == "rbf") {
    plot (D, fill=FALSE, sd=sigma);
  } else {
    plot (D, fill=FALSE);
  }
  readline ("Press any key to continue");
  G <- normalise (kernel.linear (D$X, t(D$X), type=type, sigma=sigma));
  plot (G, surface=surface, matlab=matlab, file=file);
  readline ("Press any key to continue");
}

exmp()





