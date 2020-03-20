from tools import *
import sys
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
#Script to plot the data of mpt

path = sys.argv[1]
data = sys.argv[2]


ex_mpt = mpt_data(path, [data])

if len(sys.argv) > 3:
    print(ex_mpt.guesser(sys.argv[3]))
else:
    print(ex_mpt.guesser())