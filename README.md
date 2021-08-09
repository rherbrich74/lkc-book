# Source Code for the Learning Kernel Classifiers - Theory and Algorithms Book

This repository contains the source code for all figures and a demo library of functions for classification learning. More details about the book can be found at [here](https://herbrich.me/learning-kernel-classifiers/).

## Figures

Most of the figures from the book have been created using the language [R](https://www.r-project.org/). Note that the codes are designed to output the graphics directly into a postscript file.  Moreover, the output is done into the directory `../../ps/`.

## Learning Kernel Classifiers Toolbox
In addition to the R source code of the figures, I have implemented all algorithms mentioned in Appendix D in R. The PR_LOQO routine was written by [Alexander Smola](https://alex.smola.org/) when at GMD FIRST. In order to use this library, you are required to first compile all the C source codes. On the command line, just type

```
R CMD SHLIB bayes.c
R CMD SHLIB kernels.c
R CMD SHLIB kperc.c
R CMD SHLIB pr_loqo.c
```

In order to the code in R, switch to the main directory and run `source("linear/linear.R")`. Then you can try the examples as follows:
* `exmp.bayes` and `exmp.bayes2` for Bayes point machines.
* `exmp.svm` for Support Vector machines.
* `exmp.nu.svm` for υ-support vector learning.
* `exmp.GP.class` for Gaussian process classification.
* `exmp.GP.regress` for Gaussian process regression.
* `exmp.RVM.class` for classification learning with relevance vector machines.
* `exmp.RVM.regress` for regression estimation with relevance vector machines.
* `exmp.fisher` for Fisher discriminants.
* `exmp.perc` for perceptron learning.

## Disclaimer
(C) Ralf Herbrich. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS FREE ONLY FOR NON-COMMERCIAL USE. IT MUST NOT BE MODIFIED AND DISTRIBUTED WITHOUT PRIOR PERMISSION OF THE AUTHOR. THIS SOFTWARE IS PROVIDED BY  "AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
