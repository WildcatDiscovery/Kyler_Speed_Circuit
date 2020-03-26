from tools import *
import sys
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
#Script to plot the data of mpt

path = sys.argv[1]
data = sys.argv[2]

data_edit = data.strip('\n')
sys.argv[2] = "/" + data_edit


ex_mpt = mpt_data(path, [sys.argv[2]])

if len(sys.argv) > 3:
    if len(sys.argv) > 4:
        masked_mpt = mpt_data(path, [sys.argv[2]], mask = sys.argv[4])
        print(masked_mpt.guesser())
    else:
        print(ex_mpt.guesser(sys.argv[3]))
else:
    print(ex_mpt.guesser())