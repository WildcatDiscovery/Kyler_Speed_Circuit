---
description: >-
  Some of the functions in the plot don't work if the fitting function has not
  been performed. Fitting allows us to apply the circuit equation.
---

# Fitting the Data

In a nutshell, the fitting process is an iterator, which takes in an initial set of arbitrary values, and tries different guesses out and picks the one with the lowest error. It then returns the corresponding set of coefficients. That set of coefficients can go one of two ways: it can either satisfy the threshold, which would entail that there is a combined error of less than 1e-10, or not satisfy the threshold, in which case it would be inserted back into the iterator to create a new set of coefficients.

![Basic flowchart of the fitting process.\(Assumes no masking of data\)](.gitbook/assets/image%20%2836%29.png)

This process would run until the threshold is achieved, or if this threshold is run 1000 times, the result at the end of the 1000th iteration would be it's final result. 

![Initial guesses from the introduction and the results will follow below](.gitbook/assets/image%20%2838%29.png)

```text
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 15
    # data points      = 240
    # variables        = 7
    chi-square         = 309619.885
    reduced chi-square = 1328.84071
    Akaike info crit   = 1732.99081
    Bayesian info crit = 1757.35529
[[Variables]]
    Rs:   0.40000000 +/- 39.5482383 (9887.06%) (init = 0.4)
    R:    3966.77131 +/- 99.5871306 (2.51%) (init = 3966.771)
    n:    0.65000000 +/- 0.02376492 (3.66%) (init = 0.65)
    fs:   11354.5503 +/- 729.409409 (6.42%) (init = 11354.55)
    R2:   2586062.43 +/- 399481.679 (15.45%) (init = 2586062)
    n2:   0.81097575 +/- 0.00354279 (0.44%) (init = 0.8109757)
    fs2:  0.02156687 +/- 0.00448898 (20.81%) (init = 0.02156687)
[[Correlations]] (unreported correlations are < 0.100)
    C(R2, fs2) = -0.997
    C(R, n)    = -0.862
    C(n2, fs2) =  0.851
    C(Rs, n)   =  0.814
    C(R2, n2)  = -0.810
    C(Rs, R)   = -0.782
    C(R, n2)   =  0.671
    C(R, fs2)  =  0.547
    C(R, R2)   = -0.525
    C(n, n2)   = -0.502
    C(fs, n2)  = -0.464
    C(R, fs)   = -0.444
    C(n, fs2)  = -0.413
    C(n, R2)   =  0.399
    C(fs, fs2) = -0.379
    C(fs, R2)  =  0.365
    C(Rs, n2)  = -0.353
    C(n, fs)   =  0.302
    C(Rs, fs2) = -0.291
    C(Rs, R2)  =  0.281
None
ITERATION NO:  1000
[0.4000000000005591,
 3966.7713134113706,
 0.6500000011409431,
 11354.550297315369,
 2586062.429650645,
 0.8109757489369762,
 0.021566873408073846]
```

The result above displays the initial variables at the top, and the iteration number towards the bottom and finally, the new set of coefficients. This process guesser took 1000 iterations to achieve this set of coefficients.

