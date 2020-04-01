from tools import *
import sys
import pandas as pd
#Script to display the dataframe of a single mpt file

path = sys.argv[1]
data = sys.argv[2]

data_edit = data.strip('\n')
sys.argv[2] = "/" + data_edit

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

#print(sys.argv)
print(mpt_data(path, [sys.argv[2]]).df_raw[['f', 're', 'im']])