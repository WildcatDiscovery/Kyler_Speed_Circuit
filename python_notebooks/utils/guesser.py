from tools import *
import sys
#Script to plot the data of mpt

path = sys.argv[1]
data = sys.argv[2]


ex_mpt = mpt_data(path, [data])

if len(sys.argv) > 3:
    ex_mpt.guesser(sys.argv[3])
else:
    ex_mpt.guesser()