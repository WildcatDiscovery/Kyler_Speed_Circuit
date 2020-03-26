from tools import *
import pandas as pd
import numpy as np
import statistics as stat
import sys
#Script to display the dataframe of a single mpt file

path = sys.argv[1]
data = sys.argv[2]
data_edit = data.strip('\n')
sys.argv[2] = "/" + data_edit
mask_choice = sys.argv[3]
data = sys.argv[2]


#print(sys.argv)
ex_mpt = mpt_data(path, [data])
#ex_mpt.mpt_plot()

if mask_choice == str(1):
    masker = ex_mpt.fast_mask()
    print(ex_mpt.fast_mask())
elif mask_choice == str(2):
    masker = ex_mpt.masker0()
    print(ex_mpt.masker0())
elif mask_choice == str(3):
    masker = ex_mpt.masker()
    print(ex_mpt.masker())
else:
    print("Error, not a Masking Function")
