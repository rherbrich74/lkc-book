###
### kernels.R    Kernel routines 
###
### 1999 written by Ralf Herbrich
### Technical University Berlin
###
### 2001 completed by Ralf Herbrich
### Microsoft Research Cambridge
###
### 2019 modified by Ralf Herbrich
### Amazon Development Center Germany
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

if (!is.loaded ("kinner")) {
  cat ("Loading dynamic library for Kernel computations\n")
  if (R.version$os == "Win32") {
    dyn.load ("linear/kernels.dll");
  } else {
    dyn.load ("linear/kernels.so");
  }
}

### this function setups a new kernel structure
kernel.new <- function (name, degree=1, ofs=1, sigma=1, norm=FALSE) {
  method <- switch (name,
		    "linear" = 1, 
		    "poly"   = 2,
		    "cpoly"  = 3,
		    "rbf"    = 4);
  ret <- list (type=method,
               degree=degree,
               ofs=ofs,
               sigma=sigma,
               norm=norm);
  class (ret) <- "kernel";
  return (ret);
}

### this function computes the inner products of the mapped datapoints
kernel.inner <- function (x, y, k) {
  if (class (k) != "kernel") stop ("k (3rd arg) has to be of class kernel");
  
  if (is.vector (x) && is.vector (y)) {
    if (length (x) != length (y))
      stop ("x (1st arg) and y (2nd arg) have to be of equal length");
    .C("kinner",
       as.double (x),
       as.double (y),
       as.integer (1),
       as.integer (length(x)),
       as.integer (1),
       as.integer (k$type),
       as.double (k$degree),
       as.double (k$ofs),
       as.double (k$sigma),
       as.integer (k$norm),
       ret=as.double (0))$ret;
  }
    else {
      if (is.matrix (x) && is.vector (y)) {
	if (dim (x) [2] != length (y))
	  stop ("x (1st arg) and y (2nd arg) have to be of suitable dimensions");
	matrix(.C ("kinner",
                   as.double (x),
                   as.double (y),
                   as.integer (dim (x) [1]),
                   as.integer (length (y)),
                   as.integer (1),
                   as.integer (k$type),
                   as.double (k$degree),
                   as.double (k$ofs),
                   as.double (k$sigma),
                   as.integer (k$norm),
                   ret=as.double (rep (0, dim (x) [1])))$ret, dim (x) [1]);
      }
	else {
 	  if (is.matrix (x) && is.matrix (y)) {
	    if (dim (x) [2] != dim (y) [1])
	      stop ("x (1st arg) and y (2nd arg) have to be of suitable dimensions");
	    matrix(.C ("kinner",
                       as.double (x),
                       as.double(y),
                       as.integer (dim (x) [1]),
                       as.integer (dim (x) [2]),
                       as.integer (dim (y) [2]),
                       as.integer (k$type),
                       as.double (k$degree),
                       as.double (k$ofs),
                       as.double (k$sigma),
                       as.integer (k$norm),
                       ret=as.double (rep (0, dim (x) [1] * dim (y) [2])))$ret,
                   dim (x) [1]);
	  }
	    else {
	      stop ("invalid first two arguments");
	    }
	}
    }
}

### print out a kernel
print.kernel <- function (k, ...) {
  cat ("Kernel\n-------\n\n\tType:\t", ...);
  switch (k$type,
	  cat ("linear\n", ...),
	  cat ("poly (degree: ", k$degree ,")\n", ...),
	  cat ("complete poly (degree: ", k$degree ,", ofs: ",
               k$ofs ,")\n", ...),
	  cat ("rbf  (sigma:  ", k$sigma ,")\n", ...));
  if (k$norm) {
    cat ("normalised\n");
  }
  cat ("\n", ...);
  invisible (k);
}
