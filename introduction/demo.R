###
### demo.R        produces all plots for the first chapter  
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

############################################################
## this code reads the images into memory
############################################################
read.file <- function (name="images.dat") {

  if (!file.exists (name)) {
    error ("File does not exists. Abort.");
  }

  ## read header
  header <- scan (name, n=2, quiet=TRUE);

  m <- header [1];
  N <- header [2];

  ## read the data
  D <- matrix (scan (name, skip=2, quiet=TRUE), nrow=m, byrow=TRUE);

  ## separate classes and images
  X <- D [, 1:N];
  class <- D [, (N+1)];

  return (list (X=X, class=class));
}

############################################################
## plots images from raw data
############################################################
plot.digit <- function (d, nrow=28) {
  org <- matrix (d, nrow=nrow, byrow=TRUE);
  ncol <- length(d)/nrow;
  rot <- matrix (0, nrow=nrow, ncol=ncol);
  for (i in 1:nrow) {
    rot [i, 1:ncol] <- org [nrow - i + 1, 1:ncol];
  }
  image (t (rot), zlim=c (0, 256),
         col=gray (seq (256, 0, by=-1)/256), axes=FALSE);
  box ();
}

############################################################
## creates separate postscripts for the first 100 digits from MNIST
############################################################
do.postscripts <- function () {
  D <- read.file ("mnist100.dat");
  
  cur <- rep (0, length=10);
  for (i in 1:100) {
    postscript (paste ("digit_", D$class [i], "_", cur [D$class [i] + 1],
                       ".ps", sep=""), width=1, height=1, paper="special",
                horizontal=FALSE);
    par (mai=c(0.05,0.05,0.05,0.05));
    plot.digit (D$X [i,]);
    cur [D$class [i] + 1] <- cur [D$class [i] + 1] + 1;
    dev.off ();
  }
}

############################################################
## makes function plots from a given set of points
############################################################
do.function <- function (output='SCREEN') {

  ## produce dataset
  set.seed (42);
  x <- runif (15, min=-1, max=1);
  y <- 3 * x^3 - 1.5 * x^2 - x + 3 + rnorm (15, sd=0.1);

  ## produce a uniform spacing for drawing
  x.draw <- seq (-1, 1, length=100);
  ones <- rep (1, length=length (x));
  ones.draw <- rep (1, length=100);

  ## try a linear fit
  X <- cbind (x, ones);
  w <- solve (t (X) %*% X) %*% t (X) %*% cbind (y);
  y.linear <- as.vector (cbind (x.draw, ones.draw) %*% w);

  ## try a cubic fit
  X <- cbind (x^3, x^2, x, ones);
  w <- solve (t (X) %*% X) %*% t (X) %*% cbind (y);
  y.cubic <- as.vector (cbind (x.draw^3, x.draw^2, x.draw, ones.draw) %*% w);

  ## try a cubic fit
  X <- cbind (x^10, x^9, x^8, x^7, x^6, x^5, x^4, x^3, x^2, x, ones);
  w <- solve (t (X) %*% X) %*% t (X) %*% cbind (y);
  y.tenth <- as.vector (cbind (x.draw^10, x.draw^9, x.draw^8, x.draw^7,
                               x.draw^6, x.draw^5, x.draw^4, x.draw^3,
                               x.draw^2, x.draw, ones.draw) %*% w);

  ## draw the different fits
  if (output == 'PS') {
    postscript (file="function_linear.ps");
  }
  plot (x, y, type="p", pch=4, lwd=2, cex=3, xlab="x", ylab="y",
        cex.lab=4.0, cex.axis=2.0);
  lines (x.draw, y.linear, lwd=4, lty=1);
  if (output == 'PS') {
    cat ("[function_linear.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();
  
  if (output == 'PS') {
    postscript (file="function_cubic.ps");
  }
  plot (x, y, type="p", pch=4, lwd=2, cex=3, xlab="x", ylab="y",
        cex.lab=4.0, cex.axis=2.0);
  lines (x.draw, y.cubic, lwd=4, lty=1);
  if (output == 'PS') {
    cat ("[function_cubic.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();
  
  if (output == 'PS') {
    postscript (file="function_tenth.ps");
  }
  plot (x, y, type="p", pch=4, lwd=2, cex=3, xlab="x", ylab="y",
        cex.lab=4.0, cex.axis=2.0);
  lines (x.draw, y.tenth, lwd=4, lty=1);
  if (output == 'PS') {
    cat ("[function_tenth.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();
}

############################################################
## generates a data distribution and prints some clustering
## solution as well as a estimated density
############################################################
do.unsupervised <- function (output = 'SCREEN') {

  ## produce the dataset
  set.seed (3);
  X1 <- cbind (rnorm (50, mean=+1, sd=0.1), rnorm (50, mean=0,  sd=0.1));
  X2 <- cbind (rnorm (50, mean=-1, sd=0.1), rnorm (50, mean=+1, sd=0.1));
  X3 <- cbind (rnorm (50, mean=-1, sd=0.1), rnorm (50, mean=-1, sd=0.1));
  X <- rbind (X1, X2, X3);
  
  ## generate the cluster centres
  C <- rbind (c (mean (X1 [, 1]), mean (X1 [, 2])),
              c (mean (X2 [, 1]), mean (X2 [, 2])),
              c (mean (X3 [, 1]), mean (X3 [, 2])));

  ## compute the clustering solution and plot it
  if (output == 'PS') {
    postscript (file="clustering.ps", width=5, height=5, paper="special");
  }
  x <- seq (range (X [, 1]) [1] * 1.2,
            range (X [, 1]) [2] * 1.2, length=200);
  y <- seq (range (X [, 2]) [1] * 1.2,
            range (X [, 2]) [2] * 1.2, length=200);
  d1 <- outer (x, y, function (a, b) {
    return ((a - C [1, 1])^2 + (b - C [1, 2])^2);});
  d2 <- outer (x, y, function (a, b) {
    return ((a - C [2, 1])^2 + (b - C [2, 2])^2);});
  d3 <- outer (x, y, function (a, b) {
    return ((a - C [3, 1])^2 + (b - C [3, 2])^2);});
  d.min <- pmin (d1, d2, d3);
  z <- (d.min == d1) + 2 * (d.min == d2) + 3 * (d.min == d3);
  
  plot.new ();
  plot.window (range (X [, 1]) * 1.2, range (X [, 2]) * 1.2);
  image (x, y, z, add=TRUE, col=gray (0.5 + (1:30)/100));
  points (X [, 1], X [, 2], pch=19, cex=1, col=gray (0.3));
  points (C [1, 1], C [1, 2], pch=4, lwd=2, cex=2, col=gray (1.0));
  points (C [2, 1], C [2, 2], pch=4, lwd=2, cex=2, col=gray (1.0));
  points (C [3, 1], C [3, 2], pch=4, lwd=2, cex=2, col=gray (1.0));
  box ();
  if (output == 'PS') {
    cat ("[clustering.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();

  ## compute the density estimate and plot it
  if (output == 'PS') {
    postscript (file="clustering_density.ps", width=5,
                height=5, paper="special");
  }
  x <- seq (range (X [, 1]) [1] * 1.2,
            range (X [, 1]) [2] * 1.2, length=50);
  y <- seq (range (X [, 2]) [1] * 1.2,
            range (X [, 2]) [2] * 1.2, length=50);
  d1 <- outer (x, y, function (a, b) {
    return (exp (-1/(2 * 0.07) * ((a - C [1, 1])^2 + (b - C [1, 2])^2)));});
  d2 <- outer (x, y, function (a, b) {
    return (exp (-1/(2 * 0.07) * ((a - C [2, 1])^2 + (b - C [2, 2])^2)));});
  d3 <- outer (x, y, function (a, b) {
    return (exp (-1/(2 * 0.07) * ((a - C [3, 1])^2 + (b - C [3, 2])^2)));});
  d.min <- pmin (d1, d2, d3);
  z <- (d.min == d1) + 2 * (d.min == d2) + 3 * (d.min == d3);
  
  par (mar=c(0, 0, 0, 0));
  persp (x, y, 1/3 * d1 + 1/3 * d2 + 1/3 * d3, xlab="first feature",
         ylab="second feature", zlab="density", cex=1.4,
         shade=0.75, border=gray (0.5), phi=30, theta=60);
  if (output == 'PS') {
    cat ("[clustering_density.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();
}

############################################################
## generates a data plot of the first 100 digits from the
## MNIST dataset
############################################################
do.mnist <- function (output = 'SCREEN') {
  D <- read.file ("introduction/mnist100.dat");
  
  ## draw the dataset
  if (output == 'PS') {
    postscript (file="mnist_dataset.ps", width=5, height=5,
                paper="special");
  }
  par (mfrow=c (7, 7), mar=c(0,0,0,0));
  for (i in 1:49) plot.digit (D$X [i,]);
  if (output == 'PS') {
    cat ("[mnist_dataset.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();

  ## draw the data matrix
  if (output == 'PS') {
    postscript ("mnist_dataset2.ps", width=15, height=10, paper="special");
  }
  par (mai=c (1.0,1.25,0.25,0.5));
  idx <- order (D$class [1:49]);
  image (1:(28*28), 1:49, t (256 - D$X [idx, ]), col=gray ((0:255)/255),
         xlab="features", ylab="image index", cex.lab=4.0, axes=FALSE,
         oldstyle=TRUE);
  axis (1, at=c(1, 100, 200, 300, 400, 500, 600, 700, 28*28), cex.axis=2.0);
  ticks <- 0;
  labels <- "0";
  for (j in 0:9) {
    high <- sum (D$class [1:49] <= j);
    labels <- c (labels, paste (high));
    high <- high + 0.5;
    ticks <- c (ticks, high);
    lines (c (1, 28*28), c (high, high), lwd=2);
  }
  axis (2, at=ticks, labels=labels, cex.axis=2.0);
  box ();
  if (output == 'PS') {
    cat ("[mnist_dataset2.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();

  ## draw the inner product
  if (output == 'PS') {
    postscript ("mnist_gram.ps", width=10, height=10, paper="special");
  }
  idx <- order (D$class [1:49]);
  G <- D$X [idx, ] %*% t (D$X [idx, ]);
  image (1:49, 1:49, G, col=gray ((0:255)/255),
         xlab="image index", ylab="image index",
         cex.lab=1.6, axes=FALSE, oldstyle=TRUE);
  ticks <- 0;
  labels <- "0";
  for (j in 0:9) {
    high <- sum (D$class [1:49] <= j);
    labels <- c (labels, paste (high));
    high <- high + 0.5;
    ticks <- c (ticks, high);
    lines (c (0, 50), c (high, high), lwd=2, col=gray (1.0));
    lines (c (high, high), c (0, 50), lwd=2, col=gray (1.0));
  }
  axis (1, at=ticks, labels=labels, cex.axis=1.4);
  axis (2, at=ticks, labels=labels, cex.axis=1.4);
  box ();
  if (output == 'PS') {
    cat ("[mnist_gram.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();

  ## draw the 5 closest to the 50-th digits in the first 49
  if (output == 'PS') {
    postscript ("mnist_nn.ps", width=8, height=3, paper="special");
  }
  par (mfrow=c(3,7), mar=c(0.5,0.5,0.5,0.5));
  test <- c (50, 51, 52);

  for (j in 1:length (test)) {
    dist <- double (0);
    for (i in 1:49) {
      dist <- c (dist, sqrt (sum ((D$X [i,] - D$X [test [j],])^2)));
    }
    idx <- order (dist);
    
    plot.digit (D$X [test [j],]);
    plot.new ();
    plot.digit (D$X [idx [1],]);
    plot.digit (D$X [idx [2],]);
    plot.digit (D$X [idx [3],]);
    plot.digit (D$X [idx [4],]);  
    plot.digit (D$X [idx [5],]);
  }
  if (output == 'PS') {
    cat ("[mnist_nn.ps plotted]\n");
  } else {
    readline("Press any key to continue");
  }
  dev.off ();
}

############################################################
## BOOK code
############################################################

book <- function (output = 'SCREEN') {
  ## do.postscripts ();
  do.function (output);
  do.unsupervised (output);
  do.mnist (output);
}

book()









