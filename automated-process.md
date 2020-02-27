---
description: >-
  A wrap up in all the processes used in the pages before to fully automate the
  process
---

# Automated Process

In conclusion, all the steps used above can be combined and slotted into a single function which can be used once to return a full set of optimal coefficients

```text
#Currently Testing
def full_auto(path,data):
    ex_mpt = mpt_data(path,data)
    masked_mpt = mpt_data(path,data, mask = [ex_mpt.masker()[0], ex_mpt.masker()[1]])

    Rs_guess = 40

    R_guess = 2959
    n_guess = 0.8
    fs_guess = 23023

    R2_guess = 258738
    n2_guess = 0.8
    fs2_guess = 0.2

    return masked_mpt.guesser(Rs_guess,R_guess,n_guess,fs_guess,R2_guess,n2_guess,fs2_guess)
```

```text
full_auto(r"C:\Users\cjang\Desktop\\", ['DE_40_1_30.mpt'])
```

![](.gitbook/assets/image%20%2822%29.png)

