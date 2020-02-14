#IMPORTER#
from __future__ import division
import pandas as pd
import numpy as np
from scipy.constants import codata
from pylab import *
from scipy.optimize import curve_fit
import mpmath as mp
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
#from scipy.optimize import leastsq
pd.options.mode.chained_assignment = None


def importer(path, data, cycle='off', mask=['none','none']):
        df_raw0 = []
        cycleno = []
        for j in range(len(data)):
            if data[j].find(".mpt") != -1: #file is a .mpt file
                df_raw0.append(extract_mpt(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".DTA") != -1: #file is a .dta file
                df_raw0.append(extract_dta(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".z") != -1: #file is a .z file
                df_raw0.append(extract_solar(path=path, EIS_name=data[j])) #reads all datafiles
            else:
                print('Data file(s) could not be identified')

            cycleno.append(df_raw0[j].cycle_number)
            if np.min(cycleno[j]) <= np.max(cycleno[j-1]):
                if j > 0: #corrects cycle_number except for the first data file
                    df_raw0[j].update({'cycle_number': cycleno[j]+np.max(cycleno[j-1])}) #corrects cycle number
#            else:
#                print('__init__ Error (#1)')

        #currently need to append a cycle_number coloumn to gamry files

        # adds individual dataframes into one
        if len(df_raw0) == 1:
            df_raw = df_raw0[0]
        elif len(df_raw0) == 2:
            df_raw = pd.concat([df_raw0[0], df_raw0[1]], axis=0)
        elif len(df_raw0) == 3:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2]], axis=0)
        elif len(df_raw0) == 4:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3]], axis=0)
        elif len(df_raw0) == 5:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4]], axis=0)
        elif len(df_raw0) == 6:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5]], axis=0)
        elif len(df_raw0) == 7:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6]], axis=0)
        elif len(df_raw0) == 8:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7]], axis=0)
        elif len(df_raw0) == 9:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8]], axis=0)
        elif len(df_raw0) == 10:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9]], axis=0)
        elif len(df_raw0) == 11:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9], df_raw0[10]], axis=0)
        elif len(df_raw0) == 12:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9], df_raw0[10], df_raw0[11]], axis=0)
        elif len(df_raw0) == 13:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9], df_raw0[10], df_raw0[11], df_raw0[12]], axis=0)
        elif len(df_raw0) == 14:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9], df_raw0[10], df_raw0[11]], df_raw0[12], df_raw0[13], axis=0)
        elif len(df_raw0) == 15:
            df_raw = pd.concat([df_raw0[0], df_raw0[1], df_raw0[2], df_raw0[3], df_raw0[4], df_raw0[5], df_raw0[6], df_raw0[7], df_raw0[8], df_raw0[9], df_raw0[10], df_raw0[11]], df_raw0[12], df_raw0[13], df_raw0[14], axis=0)
        else:
            print("Too many data files || 15 allowed")
        df_raw = df_raw.assign(w = 2*np.pi*df_raw.f) #creats a new coloumn with the angular frequency

        #Masking data to each cycle
        df_pre = []
        df_limited = []
        df_limited2 = []
        df = []
        if mask == ['none','none'] and cycle == 'off':
            for i in range(len(df_raw.cycle_number.unique())): #includes all data
                df.append(df_raw[df_raw.cycle_number == df_raw.cycle_number.unique()[i]])                
        elif mask == ['none','none'] and cycle != 'off':
            for i in range(len(cycle)):
                df.append(df_raw[df_raw.cycle_number == cycle[i]]) #extracting dataframe for each cycle                                
        elif mask[0] != 'none' and mask[1] == 'none' and cycle == 'off':
            df_pre = df_raw.mask(df_raw.f > mask[0])
            df_pre.dropna(how='all', inplace=True)
            for i in range(len(df_pre.cycle_number.unique())): #Appending data based on cycle number
                df.append(df_pre[df_pre.cycle_number == df_pre.cycle_number.unique()[i]])
        elif mask[0] != 'none' and mask[1] == 'none' and cycle != 'off': # or [i for i, e in enumerate(mask) if e == 'none'] == [0]
            df_limited = df_raw.mask(df_raw.f > mask[0])
            for i in range(len(cycle)):
                df.append(df_limited[df_limited.cycle_number == cycle[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle == 'off':
            df_pre = df_raw.mask(df_raw.f < mask[1])
            df_pre.dropna(how='all', inplace=True)
            for i in range(len(df_raw.cycle_number.unique())): #includes all data
                df.append(df_pre[df_pre.cycle_number == df_pre.cycle_number.unique()[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle != 'off': 
            df_limited = df_raw.mask(df_raw.f < mask[1])
            for i in range(len(cycle)):
                df.append(df_limited[df_limited.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle != 'off':
            df_limited = df_raw.mask(df_raw.f < mask[1])
            df_limited2 = df_limited.mask(df_raw.f > mask[0])
            for i in range(len(cycle)):
                df.append(df_limited[df_limited2.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle == 'off':
            df_limited = df_raw.mask(df_raw.f < mask[1])
            df_limited2 = df_limited.mask(df_raw.f > mask[0])
            for i in range(len(df_raw.cycle_number.unique())):
                df.append(df_limited[df_limited2.cycle_number == df_raw.cycle_number.unique()[i]])
        else:
            print('__init__ error (#2)')
