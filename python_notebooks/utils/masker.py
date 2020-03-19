from tools import *
import pandas as pd
import numpy as np
import statistics as stat
import sys
#Script to display the dataframe of a single mpt file

path = sys.argv[1]
data = sys.argv[2]
mask_choice = sys.argv[3]



def fast_mask(self):
        skeleton = self.df_raw.iloc[:,0:3]
        re_mid, im_mid  = np.mean(skeleton['re']), np.mean(skeleton['im'])
        a = skeleton[abs(skeleton['re']) <= re_mid * .5]
        b = skeleton[abs(skeleton['im']) <= im_mid * .5]
        c = pd.concat([a, b]).drop_duplicates()
        return [c['f'].max(), c['f'].min()]

def masker0(self):
    skeleton = self.df_raw.iloc[:,0:3]
    re_lim, im_lim  = max(skeleton['re']) * .6, max(skeleton['im'] * .6)
    a = skeleton[(skeleton['re']) <= re_lim]
    b = skeleton[(skeleton['im']) <= im_lim]
    c = pd.concat([a, b]).drop_duplicates()

    return [max(c['f']), min(c['f'])]

def masker(self, num_bins = 5):

    c = self.df_raw.iloc[:,0:3]
    res = []
    ims = []

    for i in pd.cut(c['re'], num_bins):
        res.append(i)
    for i in pd.cut(c['im'], num_bins):
        ims.append(i)
    #print('res', res)
    #print('ims', ims)
    d = c[(c['re'] >=stat.mode(res).left) & (c['re'] <= (stat.mode(res).right + (stat.mode(res).right - stat.mode(res).left)))]
    #print(stat.mode(res).left -  stat.mode(res).right)
    f = d[(d['im'] >=stat.mode(ims).left) & (d['im'] <= (stat.mode(ims).right + (stat.mode(ims).right - stat.mode(ims).left)))]
    #print(stat.mode(ims).left - stat.mode(ims).right)
    return [max(f['f']), min(f['f'])]

ex_mpt = mpt_data(path, [data])
ex_mpt.mpt_plot()

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

masked_data = mpt_data(path, [data], mask = masker)
masked_data.mpt_plot()