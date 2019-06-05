###
### demo.R          checks the tightness of the volume ratio bound
###                 by some plots
###
### 2000 written by Ralf Herbrich
### Microsoft Research Cambrdige
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

if (!is.loaded ("volratio")) {
  cat ("Loading dynamic library for volume ratio computation\n")
  if (R.version$os == "Win32") {
    dyn.load("volratio.dll");
  } else {
    dyn.load("volratio.so");    
  }
}

############################################################
## creates the bound/true value plot for a series of dimensions
############################################################

make.plots <- function (dim = c(3, 5, 7)) {
  x <- seq (0.1, 1, length=100);

  for (d in dim) {    
    exact <- vector (mode="numeric", length=length (x));
    exact <- .C ("volratio",
                 as.double (x),
                 as.integer (d),
                 as.integer (length (x)),
                 res=as.double (exact))$res;
    bd <- -d * log (x) + log (2);
    plot (range (x), range (bd, exact), type="n",
          xlab=expression (x), ylab="log volume ratio",
          cex.lab=3.5, cex.axis=2.5);
    lines (x, exact, lwd=6, lty=2);
    lines (x, bd, lwd=6, lty=1);
    legend (0.6, max (c (bd, exact))*0.95, c ("exact value", "bound"),
            lwd=4, lty=c (2, 1), cex=2.5);
    readline ("Press any key for next plot");
    dev.off ();
  }
}

############################################################
## book plots
############################################################

book <- function () {

  if (!file.exists ("../../ps/tight_d=10.ps")) {
    postscript (file="../../ps/tight_d=10.ps");
    par (mai=c (1.0,1.25,0.25,0.25));
    make.plots (dim=10);
  }

  if (!file.exists ("../../ps/tight_d=100.ps")) {
    postscript (file="../../ps/tight_d=100.ps");
    par (mai=c (1.0,1.25,0.25,0.25));
    make.plots (dim=100);
  }
}















