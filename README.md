# LegendreZeros

This repository contains MATLAB-code for a numerical scheme developed to compute the $k$-zeros of the associated Legendre function $P_{-1/2 + ik}^m(x)$, 
with $k$ and $x>1$ real and $m$ integral. Legendre functions of this form are also called Mehler functions or conical functions.

This work has been made as part of a research project at the Aalto University [Department of Mathematics and Systems Analysis](http://math.aalto.fi/en/) and 
my B.Sc. Thesis "Legendren liittofunktion nollakohdat", translating to "Zeros of the associated Legendre function".

The method uses only MATLAB files, but the open source MATLAB package Chebfun is required for proper use of the code, 
ie. the package needs to be installed and added to the MATLAB PATH. Chebfun is available for free [here](https://www.chebfun.org/) 
and the documentation can be found [here](https://www.chebfun.org/docs/).

Link to the B.Sc Thesis in the Aalto online archive [here](http://urn.fi/URN:NBN:fi:aalto-202310036125). Note that it is written in Finnish.

Instruction for use can be found in the comments of each ``.m``-file.

# Contents

* LICENSE: MIT License
* The functions ``LegendreCC``, ``LegendreIntegral`` and ``LegendreHyp`` are helper functions for evaluating the function. These are called to by the chebfun-constructor in the
* ...functions ``LegendreZeros``, ``LegendreZerosF`` and ``LegendreZero`` which are called to compute the zero(s) with specified input parameters $(m, n, x)$.
* The function ``NthAbsoluteZero`` is a simple method based on the aforementioned rootfinding methods. It computes the $n$:th $k$-zeros among all orders $m= 0, 1, 2, \dots$, with only a fixed main argument $z>1$ given as input in addition to $n$. The naming convention here is that the $k$:th _absolute_ root is $k$:th among all orders $m = 0, 1, 2, \dots$ for a fixed $z>1$. 
* The JSON-script ``functionSignatures.json`` enables input suggestions in the MATLAB editor.
