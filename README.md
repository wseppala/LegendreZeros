# LegendreZeros

This repository contains MATLAB-code for a numerical scheme developed to compute the $k$-zeros of the associated Legendre function $P_{-1/2 + ik}^m(x)$, 
with $k$ and $x>1$ real and $m$ integral. Legendre functions of this form are also called Mehler functions or conical functions.

This work has been made as part of a research project at the Aalto University [Department of Mathematics and Systems Analysis](http://math.aalto.fi/en/) and 
my B.Sc. Thesis "Legendren liittofunktion nollakohdat", translating to "Zeros of the associated Legendre function".
* NOTE: Link to the B.Sc. thesis will be added later.

The method uses only MATLAB files, but the open source MATLAB package Chebfun is required for proper use of the code, 
ie. the package needs to be installed and added to the MATLAB PATH. Chebfun is available for free [here](https://www.chebfun.org/) 
and the documentation can be found [here](https://www.chebfun.org/docs/).

Instruction for use can be found in the comments of each ``.m``-file.

# Contents

* LICENSE: MIT License
* functions ``LegendreCC``, ``LegendreIntegral`` and ``LegendreHyp`` are helper functions for evaluating the function. These are called to by the chebfun-constructor in the
* ...functions ``LegendreZeros``, ``LegendreZerosF`` and ``LegendreZero`` which are called to compute the zero(s) with specified input parameters $(m, n, x)$.
