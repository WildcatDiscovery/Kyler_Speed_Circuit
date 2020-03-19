from tools import *
import sys
#Script to plot the data of mpt

path = sys.argv[1]
data = sys.argv[2]


ex_mpt = mpt_data(path, [data])

if len(sys.argv) == 2:
    ex_mpt.mpt_plot()
elif len(sys.argv) > 3:
    if 'fit' in sys.argv:
        ex_mpt.guesser()
        ex_mpt.mpt_plot(fitting = 'on')
        if len(sys.argv) > 7:
            ex_mpt.mpt_plot(x_window = [int(sys.argv[4]), int(sys.argv[5])], y_window = [int(sys.argv[6]), int(sys.argv[7])])
    else:
        if len(sys.argv) > 6:
            ex_mpt.mpt_plot(x_window = [int(sys.argv[3]), int(sys.argv[4])], y_window = [int(sys.argv[5]), int(sys.argv[6])])
else:
    print("Error in the plot")