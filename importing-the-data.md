---
description: >-
  Currently this library works only on mpt files. You can run the import
  function on a folder full of mpt files and store the results all on a separate
  file.(TO BE IMPLEMENTED)
---

# Importing the Data

The library uses an **mpt\_data** object, which takes in the pathway to the file\(or files\), the name of the mpt file, and the mask. The mask of the mpt file is optional, that will be addressed in a later page. The example below shows a path in my data folder, where the name of the data is called _'DE\_40\_1\_30.mpt'_

```text
#Locate the data, prepare for import of the MPT file
#Import necessary packaging
from utils.tools import *

#EXAMPLE
path=r"C:\Users\cjang\Desktop\\"
data = ['DE_40_1_30.mpt']
ex_mpt = mpt_data(path,data)
```

Press **Shift + Enter** to run the cell! We have now saved that mpt's data into a mpt\_data object called **ex\_mpt**

Once imported, the data has several attributes you can access. From the picture above, we can see that our data has been saved under the name mpt\_data. We can access the information of the mpt file through viewing it as a Pandas dataframe. 

```text
ex_mpt.df_raw
```

![](.gitbook/assets/image%20%2825%29.png)

From here, we can access various information about this mpt file such as the real values column by calling 're' on the dataframe.

```text
ex_mpt.df_raw['re']
```

![We are given a list of values in the real values column, a set of float64s](.gitbook/assets/image%20%2832%29.png)

The mpt\_data object that we have created has four main functions:

1. **ex\_mpt.mpt\_plot**
2. **ex\_mpt.mpt\_fit**
3. **ex\_mpt.Lin\_KK**
4. **ex\_mpt.masker**















