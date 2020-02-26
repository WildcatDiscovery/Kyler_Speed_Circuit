---
description: >-
  Placing a limit on the frequency may help you achieve an optimal result in the
  fitting process.
---

# Masking the Data

Even on a good dataset, we can get a big error if we take into consideration every single frequency of data. Take for example, '_DE\_40\_1\_30.mpt_' which by all accounts, is a good dataset without any unusual spikes or drops in the data. Let's take a look at if we were to import all frequencies in the file.

```text
path=r"C:\Users\cjang\Desktop\\"
data = ['DE_40_1_30.mpt']
ex_mpt = mpt_data(path,data)
ex_mpt.mpt_plot()
```

![There&apos;s a lot of points on here!](.gitbook/assets/image%20%2811%29.png)

If we take the guessing iterator and run it on this dataset, the guessing iterator has to find a fitting equation that satisfies every single point on this set, or something that is close.

```text
Rs_guess = 40

R_guess = 2959
n_guess = 0.8
fs_guess = 23023

R2_guess = 258738
n2_guess = 0.8
fs2_guess = 0.2

ex_mpt.guesser(Rs_guess,R_guess,n_guess,fs_guess,R2_guess,n2_guess,fs2_guess)
```

![](.gitbook/assets/image%20%289%29.png)

![](.gitbook/assets/image%20%2821%29.png)

![The iterator reached the limit of a thousand iterations without satisfying its threshold of 1e-10.](.gitbook/assets/image%20%2818%29.png)

If we graph this 'optimal' set of coefficients, we aren't guaranteed a great fitting graph...

```text
ex_mpt.set_new_gph_dims(8,8)
ex_mpt.mpt_plot(fitting = 'on', x_window = [0,5000], y_window = [0,5000])
```

![We can do better](.gitbook/assets/image%20%2810%29.png)

What we can do is eliminate some of the frequencies from the file to make it easier to get a more accurate graph. But how do we determine which frequencies to drop from the file?

## Linear Kramer Kronig Analysis

Running the function **ex\_mpt.LinKK\(\)** will allow you to see where your residuals are fluctuating the most. From here you can determine on your own what your boundaries should be.

```text
#Will be updated
ex_mpt.Lin_KK(plot = 'w_data')
```

![It seems like our residuals stabilize from 1.75 to 6 on the x-axis. We can use this for boundaries!](.gitbook/assets/image%20%286%29.png)

```text
#Notice how the graph shows log(f) not f. We must translate back
#so instead of 1.75 and 6, we must insert 10**1.75 and 10**6
masked_mpt = mpt_data(path,data, mask = [10**6, 10**1.75])
masked_mpt.set_new_gph_dims(8,8)
masked_mpt.mpt_plot()
```

![This is much better!](.gitbook/assets/image%20%2814%29.png)

We can then fit this graph with our guessing iterator!

```text
#Same initial coefficients
Rs_guess = 40

R_guess = 2959
n_guess = 0.8
fs_guess = 23023

R2_guess = 258738
n2_guess = 0.8
fs2_guess = 0.2

masked_mpt.guesser(Rs_guess,R_guess,n_guess,fs_guess,R2_guess,n2_guess,fs2_guess)
```

![Within 7 iterations, we achieved threshold verification!](.gitbook/assets/image%20%2813%29.png)



