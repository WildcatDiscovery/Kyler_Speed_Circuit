from tools import *
import sys
#Script to plot the data of mpt

path = sys.argv[1]
data = sys.argv[2]


ex_mpt = mpt_data(path, [data])

if len(sys.argv) > 3:
    if sys.argv[3] == 'fit':
        ex_mpt.guesser()
        ex_mpt.mpt_plot(fitting = 'on')
else:
    ex_mpt.mpt_plot()