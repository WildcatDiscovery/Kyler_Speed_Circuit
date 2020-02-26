---
description: >-
  Currently this library works only on mpt files. You can run the import
  function on a folder full of mpt files and store the results all on a separate
  file.
---

# Importing the Data

The library uses an MPT\_DATA object, which takes in the pathway to the file\(or files\), the name of the mpt file, and the mask. The mask of the mpt file is optional, that will be addressed in a later page. The example below shows a path in my data folder, where the name of the data is called 'DE\_40\_1\_30.mpt'

![Make sure you put a backslash after your path!! Not doing so will cause an error!!!](.gitbook/assets/image%20%2817%29.png)

Once imported, the data has several attributes you can access. From the picture above, we can see that our data has been saved under the name mpt\_data. We can access the information of the mpt file through viewing it as a Pandas dataframe. 

![The ones we really need for our library are the first three columns: &apos;f&apos;, &apos;re&apos;, and &apos;im&apos;](.gitbook/assets/image%20%285%29.png)

The mpt\_data object that we have created has four main functions:

1. ex\_mpt.mpt\_plot
2. ex\_mpt.mpt\_fit
3. ex\_mpt.Lin\_KK
4. ex\_mpt.masker















