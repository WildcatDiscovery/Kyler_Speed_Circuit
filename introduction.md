---
description: High level view of the functionality of Kyler's Speed Circuit
---

# Introduction

Using PyEIS library, developed custom Electrochemical Impedance Spectroscopy Kristian Knudsen Developed the PyEIS library, ALL ORIGINAL PYEIS CODE BELONGS TO HIM [https://github.com/kbknudsen/PyEIS](https://github.com/kbknudsen/PyEIS) is the link to the original python library type pip install pyeis in the command prompt to install nessecary libraries

```
conda install pandas
conda install numpy
conda install scipy
conda install pylab
conda install mpmath
conda install lmfit
conda install maplotlib
conda install seaborn
```

## Function Process

Essentially, the process of how Kyler's Speed Circuit works is that it takes in an mpt file, which is then organized into a mpt object which was implemented in python, and fits a circuit equation onto it's Nyvquist Impedance graph. 



![Here we have a Nyvquist Impedance graph that has not been fitted yet.](.gitbook/assets/image%20%284%29.png)

![This graph has been fitted; the red dots overlay our initial graph](.gitbook/assets/image%20%285%29.png)

The final result yields a set of coefficients for our circuit equation. We start with a set of arbitrary initial values and iterate to the optimal set of values. We can then export these coefficients to a txt file or excel sheet.

![Here is our initial set of coefficients, which will be run into the guessing iterator](.gitbook/assets/image%20%288%29.png)

![After 15 iterations, we achieve a total error of &amp;gt;1e-10 with this set of coefficients ](.gitbook/assets/image%20%282%29.png)





