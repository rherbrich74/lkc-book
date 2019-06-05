###
### normalise.R         experiments for checking normalisation with SVMs
###
### 2000 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


source ("linear.R")

############################################################
### makes all the experiments for one dataset
############################################################
do.dataset <- function (base, kernel, degree.list, tol=1e-3) {

  ## log information
  cat ("Processing ", base, "\n======================\n\n");

  ## construct the file name of the total results file
  curves <- paste (base, '/', base, "_degree.results", sep="");
  unlink (curves);

  ##------------------------------------------------------------
  ## main loop over differen lambda values
  ##------------------------------------------------------------
  for (degree in degree.list) {

    ## log the current lambda (on screen)
    cat ("degree= ", formatC (degree, format="f"), "\n");
    
    ## prepare the error lists and the file counter
    fcnt <- 1;
    svm.error <- double (0);
    svm2.error <- double (0);
    kernel <- kernel.new ("poly", degree=degree);

    ## construct the local result filename
    log <- paste (base, '/', base, "_degree=",
                  formatC (degree, format="f"),
                  ".error", sep="");

    ## check, if the experiment is already conducted
    if (file.exists (log)) {
      
      ## if the experiments are already conducted read the files
      cat ("Reading the log file\n");
      data <- matrix (scan (log, quiet=TRUE), ncol=3, byrow=TRUE);
      svm.error <- data [,2];
      svm2.error <- data [,3];
      fcnt <- dim (data)[1];
    } else {
      cat ("Conducting experiments\n");
      
      ## unlink the error's log file
      unlink (log);

      ##------------------------------------------------------------
      ## main loop over different data files
      ##------------------------------------------------------------
      while (TRUE) {
        train.fl <- paste (base, '/', base, "_", fcnt, ".tr", sep="");
        test.fl <- paste (base, '/', base, "_", fcnt, ".ts", sep="");
        
        ## halt, if no more data files available
        if (!file.exists (train.fl) || !file.exists (test.fl)) {
          break;
        }
        ## read the datafiles
        train <- read.dataset (train.fl);
        test  <- read.dataset (test.fl);
        m <- dim (train$data) [1];
        
        ## setup a new perceptron and learn the SVM solution
        p <- perc.new (train$data, train$class, kernel);

        ## learn an SVM
        svm <- perc.svm (p, C=1e20);
        svm.err <- perc.class.error (test$data, test$class, svm);

        ## learn an SVM (normalised)
        svm2 <- perc.svm (p, norm=TRUE, C=1e20);
        svm2.err <- perc.class.error (test$data, test$class, svm2);

        ## log the error estimates in the data files
        cat ("\t", fcnt,
             "\t", formatC (svm.err, format="f"),
             "\t", formatC (svm2.err, format="f"), "\n");
        cat ("\t", fcnt,
             "\t", formatC (svm.err, format="f"),
             "\t", formatC (svm2.err, format="f"), "\n",
             file=log, append=TRUE);
        
        ## append the test.errors to the list
        svm.error <- c (svm.error, svm.err);
        svm2.error <- c (svm2.error, svm2.err);
        
        ## increment the file counter by one 
        fcnt <- fcnt + 1;
      }
      ## decrease by one (last file was not found)!
      fcnt <- fcnt - 1;
    }
    
    ## produce the result file
    result <- paste (base, '/', base, "_degree=",
                     formatC (degree, format="f"), ".results",
                     sep=""); 
    cat (base, "\n", "------------\n\n", file=result, sep="");
    cat ("degree = ", formatC (degree, format="f"), "\n\n",
         file=result, append=TRUE)
    print (kernel, file=result, append=TRUE);
    
    cat ("averaged classification error on the test sets in %\n---------------------------------------------------\n\n", file=result, append=TRUE);
    cat ("\tSVM (no normialisation): ",
         formatC (mean (svm.error), format="f"), "+-",
         formatC (sqrt (var  (svm.error) / fcnt),
                  format="f"), "\n", 
         file=result, append=TRUE);
    cat ("\tSVM (normalisation): ",
         formatC (mean (svm2.error), format="f"), "+-",
         formatC (sqrt (var  (svm2.error) / fcnt),
                  format="f"), "\n",
         file=result, append=TRUE);

    ## log the result in the global file
    cat (formatC (degree, format="f"),
         "\t", formatC (mean (svm.error), format="f"),
         "\t", formatC (sqrt (var (svm.error) / fcnt), format="f"),
         "\t", formatC (mean (svm2.error), format="f"),
         "\t", formatC (sqrt (var (svm2.error) / fcnt), format="f"),
         "\n", file=curves, append=TRUE);
  }
  return;
}


### makes all the plot outputs
do.plots <- function (base, ex=1.8) {
  
  ## log information
  cat ("Processing ", base, "\n======================\n\n");

  res.fl <- paste (base, "/", base, "_degree.results", sep="");
  if (!file.exists (res.fl)) {
    break;
  } else {
    x11 ();
    par (mfrow = c (1, 1), fig=c(0.02,1,0,1));
    
    ## read data file
    data <- matrix (scan (res.fl, quiet=TRUE), ncol=5, byrow=TRUE);
    degree <- data [,1];
    svm.gen.error.mean <- data [,2];
    svm.gen.error.std <- data [,3];
    svm2.gen.error.mean <- data [,4];
    svm2.gen.error.std <- data [,5];
    
    ## plot the result of the SVM
    plot.errorbars (degree, svm.gen.error.mean, svm.gen.error.std,
                    xlab=expression (p),
                    ylab="generalisation error",
                    cex.lab=ex, cex.axis=ex*.8,
                    ylim = c (
                      min (c (svm.gen.error.mean - svm.gen.error.std,
                              svm2.gen.error.mean - svm2.gen.error.std)),
                      max (c (svm.gen.error.mean + svm.gen.error.std,
                              svm2.gen.error.mean + svm2.gen.error.std))),
                    lwd=2, lty=1,
                    blwd=1, ticch=19, blty=1);
    plot.errorbars (degree, svm2.gen.error.mean, svm2.gen.error.std,
                    lwd=2, lty=2,
                    blwd=1, ticch=2, blty=1, add=TRUE);
    
    ps <- readline ("Type in name of postscript: ");
    dev.print (file=ps);
    
    dev.off ();
  }
  return;
}

############################################################
## M A I N   C O D E
############################################################

do.dataset ("thyroid",  degree.list=c (2, 4, 6, 8, 10, 12, 14, 16, 18, 20), tol=1e-2);
do.plots ("thyroid");

do.dataset ("sonar",  degree.list=c (2, 4, 6, 8, 10, 12, 14, 16, 18, 20), tol=1e-2);
do.plots ("sonar");









