### evidence_demo        demonstrates the idea of evidence 
###                      maximisation over 10 classes 
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.


############################################################
## EXAMPLE code
############################################################

book <- function () {
  ## plots two different distributions of classifications
  peaked <- dbinom (0:31, 31, 0.5);
  uniform <- rep (1/32, length=32);
  
  ## does the barplot
  X <- rbind (peaked, uniform);
  colnames (X) <- c ("0", "", "", "", "", "", "", "", "", "", "", "", "",
                     "", "", "", "", "", "", "", "", "", "", "", "", "",
                     "", "", "", "", "1", "1");

  postscript ("../../ps/evidence_5.ps", paper="special",
              width=10, height=6, horizontal=FALSE);
  par (cex.axis=1.2, cex.lab=1.2, cex=1.2);
  barplot (X, beside=TRUE, col=c(gray (0.3), gray (0.8)),
           legend.text=c ("simple", "uniform"));
  
  readline ("Press any key to continue");
  dev.off ();
}








