---
description: >-
  Some of the functions in the plot don't work if the fitting function has not
  been performed. Fitting allows us to apply the circuit equation.
---

# Fitting the Data

In a nutshell, the fitting process is an iterator, which takes in an initial set of arbitrary values, and tries different guesses out and picks the one with the lowest error. It then returns the corresponding set of coefficients. That set of coefficients can go one of two ways: it can either satisfy the threshold, which would entail that there is a combined error of less than 1e-10, or not satisfy the threshold, in which case it would be inserted back into the iterator to create a new set of coefficients.

![Basic flowchart of our iteration process](.gitbook/assets/image%20%289%29.png)

This process would run until the threshold is achieved, or if this threshold is run 1000 times, the result at the end of the 1000th iteration would be it's final result. 

![Initial guesses from the introduction and the results will follow below](.gitbook/assets/image%20%2825%29.png)

```text
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 9
    # data points      = 178
    # variables        = 7
    chi-square         = 206.337611
    reduced chi-square = 1.20665270
    Akaike info crit   = 40.2959698
    Bayesian info crit = 62.5684546
[[Variables]]
    Rs:   33.4894887 +/- 3.03710050 (9.07%) (init = 33.48949)
    R:    3000.12645 +/- 16.4753661 (0.55%) (init = 3000.126)
    n:    0.83870115 +/- 0.00312185 (0.37%) (init = 0.8387011)
    fs:   22738.7917 +/- 216.722198 (0.95%) (init = 22738.79)
    R2:   112840.635 +/- 20859.2424 (18.49%) (init = 112840.6)
    n2:   0.72953729 +/- 0.00645915 (0.89%) (init = 0.7295373)
    fs2:  0.94756171 +/- 0.27710293 (29.24%) (init = 0.9475617)
[[Correlations]] (unreported correlations are < 0.100)
    C(R2, fs2) = -0.999
    C(n2, fs2) =  0.935
    C(R2, n2)  = -0.917
    C(R, fs)   = -0.819
    C(R, n2)   =  0.819
    C(R, n)    = -0.748
    C(n, fs)   =  0.694
    C(R, fs2)  =  0.680
    C(R, R2)   = -0.658
    C(fs, n2)  = -0.641
    C(n, n2)   = -0.515
    C(fs, fs2) = -0.514
    C(Rs, n)   =  0.511
    C(fs, R2)  =  0.495
    C(Rs, R)   = -0.421
    C(n, fs2)  = -0.411
    C(n, R2)   =  0.396
    C(Rs, n2)  = -0.173
    C(Rs, fs)  =  0.161
    C(Rs, fs2) = -0.134
    C(Rs, R2)  =  0.128
None
ITERATION NO:  15
total error:  -2.8677504815277644e-11
[33.489488674007156,
 3000.1264502872973,
 0.8387011469573208,
 22738.79172943603,
 112840.63536435792,
 0.7295372890618985,
 0.9475617098874135]
```

The result above displays the initial variables at the top, and the iteration number towards the bottom and finally, the new set of coefficients.

