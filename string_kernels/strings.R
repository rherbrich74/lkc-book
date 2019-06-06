### strings        demonstrates the ridge problem as well as
###                the noise problem with STRING kernels
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge


############################################################
## 1) load the plotting routines from rbf_demo
## 2) loads the compiled C code
############################################################

source ("string_kernels/rbf_demo.R")

if (!is.loaded ("string_kernel")) {
  cat ("Loading dynamic library for string kernels\n")
  if (R.version$os == "Win32") {
    dyn.load ("string_kernels/strings.dll");
  } else {
    dyn.load ("string_kernels/strings.so");
  }
}

############################################################
## invokes the string program and plots the resulting Gram matrix
## ATTENTION: this program relies on the existence of the compiled
## strings.c program!
############################################################

strings <- function (collection="test", ## name of the text collection
                     type="substring",  ## type of the string kernel
                     length=5,          ## (maximal) length of substrings
                     lambda=1) {        ## (base) weight of each feature

  ## parameter check
  if (!file.exists (collection))
    stop ("Test collection does not exist.");

  if (type != "substring" &&
      type != "subsequence" &&
      type != "fullsubstring" &&
      type != "bow") 
    stop ("Kernel types must match 'substring', 'fullsubstring', 'subsequence' or 'bow'.");

  if (length < 1)
    stop ("Length must be strictly positive.");

  if (lambda < 0)
    stop ("Negative lambda make no sense");

  ## log some info
  cat ("Document library    :", collection, "\n");
  cat ("String kernel type  :", type, "\n");
  cat ("Length of substrings:", length, "\n");
  cat ("Lambda              :", lambda, "\n");
  cat ("=============================\n\n");
  
  ## count the number of documents in the directory
  ## to compute the (correct) size of the GRAM matrix
  m <- 0;
  while (file.exists (paste (collection, "/", m+1, ".txt", sep=""))) {
    m <- m + 1;
  }
  cat ("Found", m, "documents in", collection, "\n");

  ## allocate space and start routine
  G <- matrix (0, nrow=m, ncol=m);
  G <- matrix (.C ("string_kernel",
                   as.character (collection),
                   as.integer (m),
                   as.character (type),
                   as.integer (length),
                   as.double (lambda),
                   G=as.double (G))$G, nrow=m, ncol=m);
  class (G) <- "gram";

  return (G);
}

############################################################
## OLD BOOK code   !!!!!!not informative enough!!!!!
############################################################

book.old <- function () {

  ## compute the Gram matrices
  G.substring <- strings (collection="book", type="substring",
                          length=5, lambda=0.5);
  G.full <- strings (collection="book", type="fullsubstring",
                     length=5, lambda=0.5);
  G.subsequence <- strings (collection="book", type="subsequence",
                            length=5, lambda=0.5);

  ## plot the Gram matries (unnormalised)
  plot (G.substring,   matlab="p1.m", file="substring_kernel.ps");
  plot (G.full,        matlab="p2.m", file="fullstring_kernel.ps");
  plot (G.subsequence, matlab="p3.m", file="subsequence_kernel.ps");
  
  ## plot the Gram matries (normalised)
  plot (normalise (G.substring),
        matlab="p4.m", file="normalised_substring_kernel.ps");
  plot (normalise (G.full),
        matlab="p5.m", file="normalised_fullstring_kernel.ps");
  plot (normalise (G.subsequence),
        matlab="p6.m", file="normalised_subsequence_kernel.ps");
}

############################################################
## BOOK code
############################################################

book <- function (output='SCREEN') {
  
  ## plot the Gram matries (normalised)
  G.bow <- strings (collection="string_kernels/book", type="bow",
                    length=5, lambda=0.5);

  if (output == 'PS') {
    postscript (file="bow_kernel.ps");
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (normalise (G.bow), surface=FALSE);
  lines (c (11.5, 11.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (19.5, 19.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (23.5, 23.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (11.5, 11.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (19.5, 19.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (23.5, 23.5), col=gray(1.0), lwd=3);
  if (output == 'PS') {
    cat ("bow_kernel.ps created\n");
  } else {
    readline ("Press any key to continue");
  }
  dev.off ();

  G.substring <- strings (collection="string_kernels/book", type="substring",
                          length=5, lambda=0.5);
  if (output == 'PS') {
    postscript (file="substring_kernel.ps");
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (normalise (G.substring), surface=FALSE);
  lines (c (11.5, 11.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (19.5, 19.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (23.5, 23.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (11.5, 11.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (19.5, 19.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (23.5, 23.5), col=gray(1.0), lwd=3);
  if (output == 'PS') {
    cat ("substring_kernel.ps created\n");
  } else {
    readline ("Press any key to continue");
  }
  dev.off ();
  
  G.subsequence <- strings (collection="string_kernels/book", type="subsequence",
                            length=5, lambda=0.5);
  if (output == 'PS') {
    postscript (file="subsequence_kernel.ps");
  }
  par (mai=c (1.0,1.25,0.25,0.25));
  plot (normalise (G.subsequence), surface=FALSE);
  lines (c (11.5, 11.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (19.5, 19.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c (23.5, 23.5), c(0, 34), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (11.5, 11.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (19.5, 19.5), col=gray(1.0), lwd=3);
  lines (c(0, 34), c (23.5, 23.5), col=gray(1.0), lwd=3);
  if (output == 'PS') {
    cat ("subsequence_kernel.ps created\n");
  } else {
    readline ("Press any key to continue");
  }
  dev.off ();
}

book()