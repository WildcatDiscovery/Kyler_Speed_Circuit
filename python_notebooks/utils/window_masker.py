from tools import *
import pandas as pd
import numpy as np
import statistics as stat
import sys
#Script to display the dataframe of a single mpt file

path = sys.argv[1]
data = sys.argv[2]

x_min = int(sys.argv[3])
x_max = int(sys.argv[4])

y_min = int(sys.argv[5])
y_max = int(sys.argv[6])





def window_masker(self, x_window, y_window):
        adj_re = self.df_raw[(self.df_raw['re']<x_window[1]) & (self.df_raw['re']>x_window[0])]
        adj_mpt = adj_re[(adj_re['im']<y_window[1]) & (adj_re['im']>y_window[0])]
        return [max(adj_mpt['f']), min(adj_mpt['f'])]

ex_mpt = mpt_data(path, [data])
print(window_masker(ex_mpt, x_window = [x_min, x_max], y_window = [y_min, y_max]))
ex_mpt.mpt_plot(x_window = [x_min, x_max], y_window = [y_min, y_max])