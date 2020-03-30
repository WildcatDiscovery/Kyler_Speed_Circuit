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
data = sys.argv[2]

x_min = int(sys.argv[3])
x_max = int(sys.argv[4])

y_min = int(sys.argv[5])
y_max = int(sys.argv[6])


ex_mpt = mpt_data(path, [data])
#ex_mpt.mpt_plot()

if mask_choice == str(1):
    masker = ex_mpt.fast_mask()
elif mask_choice == str(2):
    masker = ex_mpt.masker0()
elif mask_choice == str(3):
    masker = ex_mpt.masker()
else:
    print("Error, not a Masking Function")
masked_mpt = mpt_data(path, [data], mask = masker)
#masked_mpt.mpt_plot()
print(masked_mpt.df_raw[['f', 're', 'im']])