###
### loss.R           depicts the different loss functions
###
### 2001 written by Ralf Herbrich
### Microsoft Research Cambridge
###
### (c) 2001 Microsoft Corporation. Reproduced with permission. All rights reserved.

book <- function () {
  x <- seq (-2, 1, length = 1000);
  heavyside <- rep (0, length=length(x));
  heavyside [x > 0] <- 1;
  hinge <- pmin (pmax (1+x, 0), 1);
  hinge2 <- pmin (pmax (1+x/2, 0), 1);
  hinge0.5 <- pmin (pmax (1+2*x, 0), 1);
  quadratic <- (pmax ((1+x), 0))^2;
  sigmoid.1 <- 2 * x + log (1 + exp (-2 * x));
  sigmoid.0.5 <- 4 * x + log (1 + exp (-4 * x));
  sigmoid.2 <- x + log (1 + exp (-x));

  ## draw HEAVYSIDE, HINGE and QUADRATIC loss
  if (!file.exists ("../../ps/loss.ps")) {
    postscript (file="../../ps/loss.ps");
    plot (c (x, 1.2), c (quadratic, 0), type="n", lwd = 1,
          xlab = expression (-yf(x)),
          ylab = "loss",
          cex.lab = 1.8, cex.axis = 1.6);
    lines (x, heavyside, lwd=4, lty=1);
    lines (x, hinge,     lwd=4, lty=2);
    lines (x, quadratic, lwd=4, lty=3);
    
    text (1.15, 0.95, expression (bold(I)[yf(x) <= 0]), cex=1.6)
    text (x [1000], hinge [1000]+0.1, "hinge loss", cex=1.6)
    text (x [1000], quadratic [1000], "quadratic loss", cex=1.6)
    
    dev.off ();
    cat ("[../../ps/loss.ps generated]\n");
  }

  ## draw several hinge losses
  if (!file.exists ("../../ps/clipped_loss.ps")) {
    postscript (file="../../ps/clipped_loss.ps");
    plot (c (x, 1.2), c (hinge0.5, 0), type="n", lwd = 1,
          xlab = expression (-yf(x)),
          ylab = "clipped loss",
          cex.lab = 1.8, cex.axis = 1.6);
    lines (x, hinge,    lwd=4, lty=1);
    lines (x, hinge2,   lwd=4, lty=2);
    lines (x, hinge0.5, lwd=4, lty=3);
    
    text (x [500], hinge [500] + 0.1, expression (paste (tau, "=1")), cex=1.6)
    text (x [500], hinge2 [500] + 0.1, expression (paste (tau, "=2")), cex=1.6)
    text (x [450], hinge0.5 [450] + 0.05, expression (paste (tau, "=0.5")), cex=1.6)
    
    dev.off ();
    cat ("[../../ps/clipped_loss.ps generated]\n");
  }

  ## draw different sigmoid's 
  if (!file.exists ("../../ps/loss_sigmoid.ps")) {
    postscript (file="../../ps/loss_sigmoid.ps");
    par (mai=c (1.0,1.25,0.25,0.25));
  
    plot (c (x, 1.2), c (sigmoid.0.5, 0), type="n", lwd = 1,
          xlab = expression (-yf(x)),
          ylab = "loss",
          cex.lab = 3.5, cex.axis = 2.5);
    lines (x, heavyside,   lwd=4, lty=1);
    lines (x, sigmoid.1,   lwd=4, lty=2);
    lines (x, sigmoid.0.5, lwd=4, lty=3);
    lines (x, sigmoid.2,   lwd=4, lty=4);
    
    text (1.05, 0.70, expression (bold(I)[yf(x) <= 0]), cex=3)
    text (x [1000], sigmoid.1 [1000]+0.1, expression (paste (beta ,"=1")),
          cex=3)
    text (x [1000], sigmoid.0.5 [1000]-0.1, expression (paste (beta ,"=0.5")),
          cex=3)
    text (x [1000], sigmoid.2 [1000]+0.1, expression (paste (beta ,"=2")),
          cex=3)
    dev.off ();
    cat ("[../../ps/loss_sigmoid.ps generated]\n");
  }
  
  ## draw the Hing-loss induced probabilities
  if (!file.exists ("../../ps/hinge_prob.ps")) {
    postscript (file="../../ps/hinge_prob.ps");
    par (mai=c (1.0,1.25,0.25,0.25));
    
    x <- seq (-2, 2, length=1000);
    prob.hinge <- exp (-pmax (1-x, 0));
    prob.hinge.inv <- exp (-pmax (1+x, 0));
    prob.sum <- prob.hinge + prob.hinge.inv;
    
    plot (c (x, 1.2), c (prob.sum+0.1, 0), type="n", lwd = 1,
          xlab = "t",
          ylab = "loss",
          cex.lab = 3.5, cex.axis = 2.5);
    lines (x, prob.hinge,     lwd=4, lty=1);
    lines (x, prob.hinge.inv, lwd=4, lty=1);
    lines (x, prob.sum,       lwd=4, lty=2);
    lines (x, rep (1, length=1000), lwd=1, lty=3);
    
    text (x [900], prob.hinge [950]-0.1,
          expression (paste (bold (P) [paste ("Y", group ("|", "T=t",""))],
              group ("(", "+1", ")"))),
          cex=3)
    
    text (x [100], prob.hinge.inv [50]-0.1,
          expression (paste (bold (P) [paste ("Y", group ("|", "T=t",""))],
              group ("(", "-1", ")"))),
          cex=3)
    
    text (x [750], prob.sum [750]+0.025,
          expression (paste (bold (P) [paste ("Y", group ("|", "T=t",""))],
              group ("(", "-1", ")"), "+",
              bold (P) [paste ("Y", group ("|", "T=t",""))],
              group ("(", "+1", ")") )),
          cex=3)
    dev.off ();
    cat ("[../../ps/hinge_prob.ps generated]\n");
  }
}

