#ASSISTING FUNCTIONS FOR THE IMPEDANCE_DATA IPYTHON NOTEBOOK
#DERIVED FROM KRISTIAN KNUDSEN'S PYEIS REPO
#HELPING TO FIT POINTS FROM A NYQUIST PLOT IN THE FORM OF A MPT FILE

#IMPORT NESSECARY LIBRARIES
#Python dependencies
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

#Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import seaborn as sns
import matplotlib.ticker as mtick
mpl.rc('mathtext', fontset='stixsans', default='regular')
mpl.rcParams.update({'axes.labelsize':22})
mpl.rc('xtick', labelsize=16) 
mpl.rc('ytick', labelsize=16)
mpl.rc('legend',fontsize=14)

from scipy.constants import codata
F = codata.physical_constants['Faraday constant'][0]
Rg = codata.physical_constants['molar gas constant'][0]

### Importing PyEIS add-ons
#from .PyEIS_Data_extraction import *
#from .PyEIS_Lin_KK import *
#from .PyEIS_Advanced_tools import *

from utils.data_extraction import *


class mpt_data:
    def __init__(self, path, data, cycle='off', mask=['none','none']):
        self.df_raw0 = []
        self.cycleno = []
        for j in range(len(data)):
            if data[j].find(".mpt") != -1: #file is a .mpt file
                self.df_raw0.append(extract_mpt(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".DTA") != -1: #file is a .dta file
                self.df_raw0.append(extract_dta(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".z") != -1: #file is a .z file
                self.df_raw0.append(extract_solar(path=path, EIS_name=data[j])) #reads all datafiles
            else:
                print('Data file(s) could not be identified')

            self.cycleno.append(self.df_raw0[j].cycle_number)
            if np.min(self.cycleno[j]) <= np.max(self.cycleno[j-1]):
                if j > 0: #corrects cycle_number except for the first data file
                    self.df_raw0[j].update({'cycle_number': self.cycleno[j]+np.max(self.cycleno[j-1])}) #corrects cycle number
#            else:
#                print('__init__ Error (#1)')

        #currently need to append a cycle_number coloumn to gamry files

        # adds individual dataframes into one
        self.df_raw = [i for i in self.df_raw0][0]

        #Masking data to each cycle
        self.df_pre = []
        self.df_limited = []
        self.df_limited2 = []
        self.df = []
        if mask == ['none','none'] and cycle == 'off':
            for i in range(len(self.df_raw.cycle_number.unique())): #includes all data
                self.df.append(self.df_raw[self.df_raw.cycle_number == self.df_raw.cycle_number.unique()[i]])                
        elif mask == ['none','none'] and cycle != 'off':
            for i in range(len(cycle)):
                self.df.append(self.df_raw[self.df_raw.cycle_number == cycle[i]]) #extracting dataframe for each cycle                                
        elif mask[0] != 'none' and mask[1] == 'none' and cycle == 'off':
            self.df_pre = self.df_raw.mask(self.df_raw.f > mask[0])
            self.df_pre.dropna(how='all', inplace=True)
            for i in range(len(self.df_pre.cycle_number.unique())): #Appending data based on cycle number
                self.df.append(self.df_pre[self.df_pre.cycle_number == self.df_pre.cycle_number.unique()[i]])
        elif mask[0] != 'none' and mask[1] == 'none' and cycle != 'off': # or [i for i, e in enumerate(mask) if e == 'none'] == [0]
            self.df_limited = self.df_raw.mask(self.df_raw.f > mask[0])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited.cycle_number == cycle[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle == 'off':
            self.df_pre = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_pre.dropna(how='all', inplace=True)
            for i in range(len(self.df_raw.cycle_number.unique())): #includes all data
                self.df.append(self.df_pre[self.df_pre.cycle_number == self.df_pre.cycle_number.unique()[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle != 'off': 
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle != 'off':
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_limited2 = self.df_limited.mask(self.df_raw.f > mask[0])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited2.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle == 'off':
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_limited2 = self.df_limited.mask(self.df_raw.f > mask[0])
            for i in range(len(self.df_raw.cycle_number.unique())):
                self.df.append(self.df_limited[self.df_limited2.cycle_number == self.df_raw.cycle_number.unique()[i]])
        else:
            print('__init__ error (#2)')

    def mpt_plot(self, bode='off', fitting='off', rr='off', nyq_xlim='none', nyq_ylim='none', legend='on', savefig='none'):
        if bode=='off':
            fig = figure(dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(111, aspect='equal')

        elif bode=='on' and rr=='off' or bode=='log' and rr=='off' or bode=='re' and rr=='off' or bode=='log_re' and rr=='off' or bode=='im' and rr=='off' or bode=='log_im' and rr=='off' or bode=='log' and rr=='off':
            fig = figure(figsize=(6, 5), dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(211, aspect='equal')
            ax1 = fig.add_subplot(212)

        elif bode=='on' and rr=='on' or bode=='log' and rr=='on' or bode=='re' and rr=='on' or bode=='log_re' and rr=='on' or bode=='im' and rr=='on' or bode=='log_im' and rr=='on' or bode=='log' and rr=='on':
            fig = figure(figsize=(6, 8), dpi=120, facecolor='w', edgecolor='k')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(311, aspect='equal')
            ax1 = fig.add_subplot(312)
            ax2 = fig.add_subplot(313)

        ### Colors
        colors = sns.color_palette("colorblind", n_colors=len(self.df))
        colors_real = sns.color_palette("Blues", n_colors=len(self.df)+2)
        colors_imag = sns.color_palette("Oranges", n_colors=len(self.df)+2)

        ### Label functions
        self.label_re_1 = []
        self.label_im_1 = []
        self.label_cycleno = []
        if legend == 'on':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' (#"+str(i+1)+")")
                self.label_im_1.append("Z'' (#"+str(i+1)+")")
                self.label_cycleno.append('#'+str(i+1))
        elif legend == 'potential':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_im_1.append("Z'' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_cycleno.append(str(np.round(np.average(self.df[i].E_avg), 2))+' V')



        ### Nyquist Plot
        for i in range(len(self.df)):
            ax.plot(self.df[i].re, self.df[i].im, marker='o', ms=4, lw=2, color=colors[i], ls='-', label=self.label_cycleno[i])
            if fitting == 'on':
                ax.plot(self.circuit_fit[i].real, -self.circuit_fit[i].imag, lw=0, marker='o', ms=8, mec='r', mew=1, mfc='none', label='')

        ### Bode Plot
        if bode=='on':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i], marker='D', ms=3, lw=2.25, ls='-', label=self.label_re_1[i])
                ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i], marker='s', ms=3, lw=2.25, ls='-', label=self.label_im_1[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), self.circuit_fit[i].real, lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='')
                    ax1.plot(np.log10(self.df[i].f), -self.circuit_fit[i].imag, lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("Z', -Z'' [$\Omega$]")
                if legend == 'on' or legend == 'potential': 
                    ax1.legend(loc='best', fontsize=10, frameon=False)
            
        elif bode == 're':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i], marker='D', ms=3, lw=2.25, ls='-', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), self.circuit_fit[i].real, lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("Z' [$\Omega$]")
                if legend == 'on' or legend =='potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)

        elif bode == 'log_re':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i], marker='D', ms=3, lw=2.25, ls='-', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.circuit_fit[i].real), lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(Z') [$\Omega$]")
                if legend == 'on' or legend == 'potential': 
                    ax1.legend(loc='best', fontsize=10, frameon=False)

        elif bode=='im':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i], marker='s', ms=3, lw=2.25, ls='-', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), -self.circuit_fit[i].imag, lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none', label='')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("-Z'' [$\Omega$]")
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)

        elif bode=='log_im':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i], marker='s', ms=3, lw=2.25, ls='-', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(-self.circuit_fit[i].imag), lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none', label='')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(-Z'') [$\Omega$]")
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)

        elif bode == 'log':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i], marker='D', ms=3, lw=2.25,  ls='-', label=self.label_re_1[i])
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i], marker='s', ms=3, lw=2.25,  ls='-', label=self.label_im_1[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.circuit_fit[i].real), lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='')
                    ax1.plot(np.log10(self.df[i].f), np.log10(-self.circuit_fit[i].imag), lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)

    def mpt_fit(self, params, circuit, weight_func='modulus', nan_policy='raise'):
        '''
        EIS_fit() fits experimental data to an equivalent circuit model using complex non-linear least-squares (CNLS) fitting procedure and allows for batch fitting.
        
        Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
        
        Inputs
        ------------
        - circuit:
          Choose an equivalent circuits and defined circuit as a string. The following circuits are avaliable.
            - RC
            - RQ
            - R-RQ
            - R-RQ-RQ
            - R-Q
            - R-RQ-Q
            - R-(Q(RW))
            - C-RC-C
            - Q-RQ-Q
            - RC-RC-ZD
            - R-TLsQ
            - R-RQ-TLsQ
            - R-TLs
            - R-RQ-TLs
            - R-TLQ
            - R-RQ-TLQ
            - R-TL
            - R-RQ-TL
            - R-TL1Dsolid (reactive interface with 1D solid-state diffusion)
            - R-RQ-TL1Dsolid

        - weight_func
          The weight function to which the CNLS fitting is performed
            - modulus (default)
            - unity
            - proportional
        
        - nan_policy
        How to handle Nan or missing values in dataset
            - ‘raise’ = raise a value error (default)
            - ‘propagate’ = do nothing
            - ‘omit’ = drops missing data
        
        Returns
        ------------
        Returns the fitted impedance spectra(s) but also the fitted parameters that were used in the initial guesses. To call these use e.g. self.fit_Rs
        '''
        self.Fit = []
        self.circuit_fit = []
        self.fit_E = []
        for i in range(len(self.df)):
            self.Fit.append(minimize(leastsq_errorfunc, params, method='leastsq', args=(self.df[i].w.values, self.df[i].re.values, self.df[i].im.values, circuit, weight_func), nan_policy=nan_policy, maxfev=9999990))
            print(report_fit(self.Fit[i]))
            
            self.fit_E.append(np.average(self.df[i].E_avg))
            
        if circuit == 'C':
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(elem_C(w=self.df[i].w, C=self.Fit[i].params.get('C').value))
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'Q':
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(elem_Q(w=self.df[i].w, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'R-C':
            self.fit_Rs = []
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'R-Q':
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'RC':
            self.fit_R = []
            self.fit_C = []
            self.fit_fs = []
            for i in range(len(self.df)):
                if "'C'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RC(w=self.df[i].w, C=self.Fit[i].params.get('C').value, R=self.Fit[i].params.get('R').value, fs='none'))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                elif "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RC(w=self.df[i].w, C='none', R=self.Fit[i].params.get('R').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_fs.append(self.Fit[i].params.get('R').value)
        elif circuit == 'RQ':
            self.fit_R = []
            self.fit_n = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RQ(w=self.df[i].w, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RQ(w=self.df[i].w, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
        elif circuit == 'R-RQ':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
        elif circuit == 'R-RQ-RQ':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_fs = []
            self.fit_fs2 = []
            self.fit_Q = []
            self.fit_Q2 = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value, R2=self.Fit[i].params.get('R2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                elif "'Q'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none', R2=self.Fit[i].params.get('R2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                elif "'fs'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, fs2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                elif "'Q'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none', R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, fs2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
        elif circuit == 'R-RC-C':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_C1 = []
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsRCC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, C1=self.Fit[i].params.get('C1').value, C=self.Fit[i].params.get('C').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R1.append(self.Fit[i].params.get('R1').value)
                self.fit_C1.append(self.Fit[i].params.get('C1').value)
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'R-RC-Q':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_C1 = []
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsRCQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, C1=self.Fit[i].params.get('C1').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R1.append(self.Fit[i].params.get('R1').value)
                self.fit_C1.append(self.Fit[i].params.get('C1').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'R-RQ-Q':
            self.fit_Rs = []
            self.fit_n = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R1=self.Fit[i].params.get('R1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, fs1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
        elif circuit == 'R-RQ-C':
            self.fit_Rs = []
            self.fit_C = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q1 = []
            self.fit_fs1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value, R1=self.Fit[i].params.get('R1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, fs1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
        elif circuit == 'R-(Q(RW))':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_sigma = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Randles_simplified(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, fs='none', n=self.Fit[i].params.get('n').value, sigma=self.Fit[i].params.get('sigma').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_sigma.append(self.Fit[i].params.get('sigma').value)
                elif "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Randles_simplified(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', fs=self.Fit[i].params.get('fs').value, n=self.Fit[i].params.get('n').value, sigma=self.Fit[i].params.get('sigma').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_sigma.append(self.Fit[i].params.get('sigma').value)
        elif circuit == 'R-TLsQ':
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-RQ-TLsQ':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, n1=self.Fit[i].params.get('n1').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Q1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', n1=self.Fit[i].params.get('n1').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Q1=self.Fit[i].params.get('Q1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-TLs':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-RQ-TLs':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs1 = []
            self.fit_fs2 = []
            self.fit_Q1 = []
            self.fit_Q2 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value, Q1='none', Q2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value, Q1=self.Fit[i].params.get('Q1').value, Q2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2='none', Q1='none', Q2=self.Fit[i].params.get('Q2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2='none', Q1=self.Fit[i].params.get('Q1').value, Q2=self.Fit[i].params.get('Q2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-TLQ':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                self.fit_L.append(self.Fit[i].params.get('L').value)            
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)            
                self.fit_Q.append(self.Fit[i].params.get('Q').value)            
                self.fit_n.append(self.Fit[i].params.get('n').value)            
                self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TLQ':
            self.fit_Rs = []
            self.fit_L = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none'))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
        elif circuit == 'R-TL':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_R = []
            self.fit_fs = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, fs=self.Fit[i].params.get('fs').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, Q='none'))                
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                
                    self.fit_R.append(self.Fit[i].params.get('R').value)                
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)                
                    self.fit_n.append(self.Fit[i].params.get('n').value)                
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TL':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_Q1 = []
            self.fit_Q2 = []
            self.fit_fs1 = []
            self.fit_fs2 = []
            for i in range(len(self.df)):
                if "'Q1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2='none', Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                    
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)                    
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)                    
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)                    
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)                    
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2=self.Fit[i].params.get('fs2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2=self.Fit[i].params.get('fs2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2='none', Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-TL1Dsolid':
            self.fit_L = []
            self.fit_radius = []
            self.fit_D = []
            self.fit_Rs = []
            self.fit_R = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_R_w = []
            self.fit_n_w = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_RsTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                self.fit_L.append(self.Fit[i].params.get('L').value)
                self.fit_radius.append(self.Fit[i].params.get('radius').value)
                self.fit_D.append(self.Fit[i].params.get('D').value)            
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R.append(self.Fit[i].params.get('R').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
                self.fit_R_w.append(self.Fit[i].params.get('R_w').value)
                self.fit_n_w.append(self.Fit[i].params.get('n_w').value)
                self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TL1Dsolid':
            self.fit_L = []
            self.fit_radius = []
            self.fit_D = []
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_Q2 = []
            self.fit_n2 = []
            self.fit_R_w = []
            self.fit_n_w = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, Q1='none', fs1=self.Fit[i].params.get('fs1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                    
                    self.fit_radius.append(self.Fit[i].params.get('radius').value)                    
                    self.fit_D.append(self.Fit[i].params.get('D').value)                                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)                    
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)                    
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)                    
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)                    
                    self.fit_R_w.append(self.Fit[i].params.get('R_w').value)                    
                    self.fit_n_w.append(self.Fit[i].params.get('n_w').value)                    
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_RsRQTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, fs1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_radius.append(self.Fit[i].params.get('radius').value)
                    self.fit_D.append(self.Fit[i].params.get('D').value)            
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_R_w.append(self.Fit[i].params.get('R_w').value)
                    self.fit_n_w.append(self.Fit[i].params.get('n_w').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'C-RC-C':
            self.fit_Ce = []
            self.fit_Rb = []
            self.fit_fsb = []
            self.fit_Cb = []
            for i in range(len(self.df)):
                if "'fsb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_C_RC_C(w=self.df[i].w, Ce=self.Fit[i].params.get('Ce').value, Cb='none', Rb=self.Fit[i].params.get('Rb').value, fsb=self.Fit[i].params.get('fsb').value))                    
                    self.fit_Ce.append(self.Fit[i].params.get('Ce').value)                    
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)
                    self.fit_fsb.append(self.Fit[i].params.get('fsb').value)
                elif "'Cb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_C_RC_C(w=self.df[i].w, Ce=self.Fit[i].params.get('Ce').value, Cb=self.Fit[i].params.get('Cb').value, Rb=self.Fit[i].params.get('Rb').value, fsb='none'))
                    self.fit_Ce.append(self.Fit[i].params.get('Ce').value)
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)
                    self.fit_Cb.append(self.Fit[i].params.get('Cb').value)
        elif circuit == 'Q-RQ-Q':
            self.fit_Qe = []
            self.fit_ne = []
            self.fit_Rb = []
            self.fit_nb = []
            self.fit_fsb = []
            self.fit_Qb = []
            for i in range(len(self.df)):
                if "'fsb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Q_RQ_Q(w=self.df[i].w, Qe=self.Fit[i].params.get('Qe').value, ne=self.Fit[i].params.get('ne').value, Qb='none', Rb=self.Fit[i].params.get('Rb').value, fsb=self.Fit[i].params.get('fsb').value, nb=self.Fit[i].params.get('nb').value))
                    self.fit_Qe.append(self.Fit[i].params.get('Qe').value)                    
                    self.fit_ne.append(self.Fit[i].params.get('ne').value)                    
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)                    
                    self.fit_fsb.append(self.Fit[i].params.get('fsb').value)
                    self.fit_nb.append(self.Fit[i].params.get('nb').value)
                elif "'Qb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Q_RQ_Q(w=self.df[i].w, Qe=self.Fit[i].params.get('Qe').value, ne=self.Fit[i].params.get('ne').value, Qb=self.Fit[i].params.get('Qb').value, Rb=self.Fit[i].params.get('Rb').value, fsb='none', nb=self.Fit[i].params.get('nb').value))
                    self.fit_Qe.append(self.Fit[i].params.get('Qe').value)
                    self.fit_ne.append(self.Fit[i].params.get('ne').value)
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)                    
                    self.fit_Qb.append(self.Fit[i].params.get('Qb').value)
                    self.fit_nb.append(self.Fit[i].params.get('nb').value)
        else:
            print('Circuit was not properly defined, see details described in definition')

    def Lin_KK(self, num_RC='auto', legend='on', plot='residuals', bode='off', nyq_xlim='none', nyq_ylim='none', weight_func='Boukamp', savefig='none'):
        '''
        Plots the Linear Kramers-Kronig (KK) Validity Test
        The script is based on Boukamp and Schōnleber et al.'s papers for fitting the resistances of multiple -(RC)- circuits
        to the data. A data quality analysis can hereby be made on the basis of the relative residuals

        Ref.:
            - Schōnleber, M. et al. Electrochimica Acta 131 (2014) 20-27
            - Boukamp, B.A. J. Electrochem. Soc., 142, 6, 1885-1894 
        
        The function performs the KK analysis and as default the relative residuals in each subplot        
    
        Note, that weigh_func should be equal to 'Boukamp'.
        
        Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)
        
        Optional Inputs
        -----------------
        - num_RC:
            - 'auto' applies an automatic algorithm developed by Schōnleber, M. et al. Electrochimica Acta 131 (2014) 20-27
            that ensures no under- or over-fitting occurs
            - can be hardwired by inserting any number (RC-elements/decade)

        - plot: 
            - 'residuals' = plots the relative residuals in subplots correspoding to the cycle numbers picked
            - 'w_data' = plots the relative residuals with the experimental data, in Nyquist and bode plot if desired, see 'bode =' in description
        
        - nyq_xlim/nyq_xlim: Change the x/y-axis limits on nyquist plot, if not equal to 'none' state [min,max] value
        
        - legend:
            - 'on' = displays cycle number
            - 'potential' = displays average potential which the spectra was measured at
            - 'off' = off

        bode = Plots Bode Plot - options:
            'on' = re, im vs. log(freq)
            'log' = log(re, im) vs. log(freq)
            
            're' = re vs. log(freq)
            'log_re' = log(re) vs. log(freq)
            
            'im' = im vs. log(freq)
            'log_im' = log(im) vs. log(freq)
        '''
        if num_RC == 'auto':
            print('cycle || No. RC-elements ||   u')
            self.decade = []
            self.Rparam = []
            self.t_const = []
            self.Lin_KK_Fit = []
            self.R_names = []
            self.KK_R0 = []
            self.KK_R = []
            self.number_RC = []
            self.number_RC_sort = []
    
            self.KK_u = []
            self.KK_Rgreater = []
            self.KK_Rminor = []
            M = 2
            for i in range(len(self.df)):
                self.decade.append(np.log10(np.max(self.df[i].f))-np.log10(np.min(self.df[i].f))) #determine the number of RC circuits based on the number of decades measured and num_RC
                self.number_RC.append(M)
                self.number_RC_sort.append(M) #needed for self.KK_R
                self.Rparam.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[0]) #Creates intial guesses for R's
                self.t_const.append(KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC[i]))) #Creates time constants values for self.number_RC -(RC)- circuits
                
                self.Lin_KK_Fit.append(minimize(KK_errorfunc, self.Rparam[i], method='leastsq', args=(self.df[i].w.values, self.df[i].re.values, self.df[i].im.values, self.number_RC[i], weight_func, self.t_const[i]) )) #maxfev=99
                self.R_names.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[1]) #creates R names
                for j in range(len(self.R_names[i])):
                    self.KK_R0.append(self.Lin_KK_Fit[i].params.get(self.R_names[i][j]).value)
            self.number_RC_sort.insert(0,0) #needed for self.KK_R
            for i in range(len(self.df)):
                self.KK_R.append(self.KK_R0[int(np.cumsum(self.number_RC_sort)[i]):int(np.cumsum(self.number_RC_sort)[i+1])]) #assigns resistances from each spectra to their respective df
                self.KK_Rgreater.append(np.where(np.array(self.KK_R)[i] >= 0, np.array(self.KK_R)[i], 0) )
                self.KK_Rminor.append(np.where(np.array(self.KK_R)[i] < 0, np.array(self.KK_R)[i], 0) )
                self.KK_u.append(1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i]))))
            
            for i in range(len(self.df)):
                while self.KK_u[i] <= 0.75 or self.KK_u[i] >= 0.88:
                    self.number_RC_sort0 = []
                    self.KK_R_lim = []
                    self.number_RC[i] = self.number_RC[i] + 1
                    self.number_RC_sort0.append(self.number_RC)
                    self.number_RC_sort = np.insert(self.number_RC_sort0, 0,0)
                    self.Rparam[i] = KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[0] #Creates intial guesses for R's
                    self.t_const[i] = KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC[i])) #Creates time constants values for self.number_RC -(RC)- circuits
                    self.Lin_KK_Fit[i] = minimize(KK_errorfunc, self.Rparam[i], method='leastsq', args=(self.df[i].w.values, self.df[i].re.values, self.df[i].im.values, self.number_RC[i], weight_func, self.t_const[i]) ) #maxfev=99
                    self.R_names[i] = KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[1] #creates R names
                    self.KK_R0 = np.delete(np.array(self.KK_R0), np.s_[0:len(self.KK_R0)])
                    self.KK_R0 = []
                    for q in range(len(self.df)):
                        for j in range(len(self.R_names[q])):
                            self.KK_R0.append(self.Lin_KK_Fit[q].params.get(self.R_names[q][j]).value)
                    self.KK_R_lim = np.cumsum(self.number_RC_sort) #used for KK_R[i]
    
                    self.KK_R[i] = self.KK_R0[self.KK_R_lim[i]:self.KK_R_lim[i+1]] #assigns resistances from each spectra to their respective df
                    self.KK_Rgreater[i] = np.where(np.array(self.KK_R[i]) >= 0, np.array(self.KK_R[i]), 0)
                    self.KK_Rminor[i] = np.where(np.array(self.KK_R[i]) < 0, np.array(self.KK_R[i]), 0)
                    self.KK_u[i] = 1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i])))
                else:
                    print('['+str(i+1)+']'+'            '+str(self.number_RC[i]),'           '+str(np.round(self.KK_u[i],2)))

        elif num_RC != 'auto': #hardwired number of RC-elements/decade
            print('cycle ||   u')
            self.decade = []
            self.number_RC0 = []
            self.number_RC = []
            self.Rparam = []
            self.t_const = []
            self.Lin_KK_Fit = []
            self.R_names = []
            self.KK_R0 = []
            self.KK_R = []
            for i in range(len(self.df)):
                self.decade.append(np.log10(np.max(self.df[i].f))-np.log10(np.min(self.df[i].f))) #determine the number of RC circuits based on the number of decades measured and num_RC
                self.number_RC0.append(np.round(num_RC * self.decade[i]))
                self.number_RC.append(np.round(num_RC * self.decade[i])) #Creats the the number of -(RC)- circuits
                self.Rparam.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC0[i]))[0]) #Creates intial guesses for R's
                self.t_const.append(KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC0[i]))) #Creates time constants values for self.number_RC -(RC)- circuits
                self.Lin_KK_Fit.append(minimize(KK_errorfunc, self.Rparam[i], method='leastsq', args=(self.df[i].w.values, self.df[i].re.values, self.df[i].im.values, self.number_RC0[i], weight_func, self.t_const[i]) )) #maxfev=99
                self.R_names.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC0[i]))[1]) #creates R names            
                for j in range(len(self.R_names[i])):
                    self.KK_R0.append(self.Lin_KK_Fit[i].params.get(self.R_names[i][j]).value)
            self.number_RC0.insert(0,0)
    
    #        print(report_fit(self.Lin_KK_Fit[i])) # prints fitting report
    
            self.KK_circuit_fit = []
            self.KK_rr_re = []
            self.KK_rr_im = []
            self.KK_Rgreater = []
            self.KK_Rminor = []
            self.KK_u = []
            for i in range(len(self.df)):
                self.KK_R.append(self.KK_R0[int(np.cumsum(self.number_RC0)[i]):int(np.cumsum(self.number_RC0)[i+1])]) #assigns resistances from each spectra to their respective df
                self.KK_Rx = np.array(self.KK_R)
                self.KK_Rgreater.append(np.where(self.KK_Rx[i] >= 0, self.KK_Rx[i], 0) )
                self.KK_Rminor.append(np.where(self.KK_Rx[i] < 0, self.KK_Rx[i], 0) )
                self.KK_u.append(1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i])))) #currently gives incorrect values
                print('['+str(i+1)+']'+'       '+str(np.round(self.KK_u[i],2)))
        else:
            print('num_RC incorrectly defined')

        self.KK_circuit_fit = []
        self.KK_rr_re = []
        self.KK_rr_im = []
        for i in range(len(self.df)):
            if int(self.number_RC[i]) == 2:
                self.KK_circuit_fit.append(KK_RC2(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 3:
                self.KK_circuit_fit.append(KK_RC3(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 4:
                self.KK_circuit_fit.append(KK_RC4(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 5:
                self.KK_circuit_fit.append(KK_RC5(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 6:
                self.KK_circuit_fit.append(KK_RC6(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 7:
                self.KK_circuit_fit.append(KK_RC7(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 8:
                self.KK_circuit_fit.append(KK_RC8(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 9:
                self.KK_circuit_fit.append(KK_RC9(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 10:
                self.KK_circuit_fit.append(KK_RC10(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 11:
                self.KK_circuit_fit.append(KK_RC11(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 12:
                self.KK_circuit_fit.append(KK_RC12(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 13:
                self.KK_circuit_fit.append(KK_RC13(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 14:
                self.KK_circuit_fit.append(KK_RC14(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 15:
                self.KK_circuit_fit.append(KK_RC15(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 16:
                self.KK_circuit_fit.append(KK_RC16(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 17:
                self.KK_circuit_fit.append(KK_RC17(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 18:
                self.KK_circuit_fit.append(KK_RC18(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 19:
                self.KK_circuit_fit.append(KK_RC19(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 20:
                self.KK_circuit_fit.append(KK_RC20(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 21:
                self.KK_circuit_fit.append(KK_RC21(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 22:
                self.KK_circuit_fit.append(KK_RC22(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 23:
                self.KK_circuit_fit.append(KK_RC23(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 24:
                self.KK_circuit_fit.append(KK_RC24(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 25:
                self.KK_circuit_fit.append(KK_RC25(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 26:
                self.KK_circuit_fit.append(KK_RC26(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 27:
                self.KK_circuit_fit.append(KK_RC27(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 28:
                self.KK_circuit_fit.append(KK_RC28(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 29:
                self.KK_circuit_fit.append(KK_RC29(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 30:
                self.KK_circuit_fit.append(KK_RC30(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 31:
                self.KK_circuit_fit.append(KK_RC31(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 32:
                self.KK_circuit_fit.append(KK_RC32(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 33:
                self.KK_circuit_fit.append(KK_RC33(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 34:
                self.KK_circuit_fit.append(KK_RC34(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 35:
                self.KK_circuit_fit.append(KK_RC35(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 36:
                self.KK_circuit_fit.append(KK_RC36(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 37:
                self.KK_circuit_fit.append(KK_RC37(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 38:
                self.KK_circuit_fit.append(KK_RC38(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 39:
                self.KK_circuit_fit.append(KK_RC39(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 40:
                self.KK_circuit_fit.append(KK_RC40(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 41:
                self.KK_circuit_fit.append(KK_RC41(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 42:
                self.KK_circuit_fit.append(KK_RC42(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 43:
                self.KK_circuit_fit.append(KK_RC43(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 44:
                self.KK_circuit_fit.append(KK_RC44(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 45:
                self.KK_circuit_fit.append(KK_RC45(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 46:
                self.KK_circuit_fit.append(KK_RC46(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 47:
                self.KK_circuit_fit.append(KK_RC47(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 48:
                self.KK_circuit_fit.append(KK_RC48(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 49:
                self.KK_circuit_fit.append(KK_RC49(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 50:
                self.KK_circuit_fit.append(KK_RC50(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 51:
                self.KK_circuit_fit.append(KK_RC51(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 52:
                self.KK_circuit_fit.append(KK_RC52(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 53:
                self.KK_circuit_fit.append(KK_RC53(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 54:
                self.KK_circuit_fit.append(KK_RC54(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 55:
                self.KK_circuit_fit.append(KK_RC55(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 56:
                self.KK_circuit_fit.append(KK_RC56(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 57:
                self.KK_circuit_fit.append(KK_RC57(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 58:
                self.KK_circuit_fit.append(KK_RC58(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 59:
                self.KK_circuit_fit.append(KK_RC59(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 60:
                self.KK_circuit_fit.append(KK_RC60(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 61:
                self.KK_circuit_fit.append(KK_RC61(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 62:
                self.KK_circuit_fit.append(KK_RC62(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 63:
                self.KK_circuit_fit.append(KK_RC63(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 64:
                self.KK_circuit_fit.append(KK_RC64(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 65:
                self.KK_circuit_fit.append(KK_RC65(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 66:
                self.KK_circuit_fit.append(KK_RC66(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 67:
                self.KK_circuit_fit.append(KK_RC67(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 68:
                self.KK_circuit_fit.append(KK_RC68(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 69:
                self.KK_circuit_fit.append(KK_RC69(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 70:
                self.KK_circuit_fit.append(KK_RC70(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 71:
                self.KK_circuit_fit.append(KK_RC71(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 72:
                self.KK_circuit_fit.append(KK_RC72(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 73:
                self.KK_circuit_fit.append(KK_RC73(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 74:
                self.KK_circuit_fit.append(KK_RC74(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 75:
                self.KK_circuit_fit.append(KK_RC75(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 76:
                self.KK_circuit_fit.append(KK_RC76(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 77:
                self.KK_circuit_fit.append(KK_RC77(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 78:
                self.KK_circuit_fit.append(KK_RC78(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 79:
                self.KK_circuit_fit.append(KK_RC79(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 80:
                self.KK_circuit_fit.append(KK_RC80(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            else:
                print('RC simulation circuit not defined')
                print('   Number of RC = ', self.number_RC)
            self.KK_rr_re.append(residual_real(re=self.df[i].re, fit_re=self.KK_circuit_fit[i].real, fit_im=-self.KK_circuit_fit[i].imag)) #relative residuals for the real part
            self.KK_rr_im.append(residual_imag(im=self.df[i].im, fit_re=self.KK_circuit_fit[i].real, fit_im=-self.KK_circuit_fit[i].imag)) #relative residuals for the imag part

        ### Plotting Linear_kk results
        ##
        #
        ### Label functions
        self.label_re_1 = []
        self.label_im_1 = []
        self.label_cycleno = []
        if legend == 'on':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' (#"+str(i+1)+")")
                self.label_im_1.append("Z'' (#"+str(i+1)+")")
                self.label_cycleno.append('#'+str(i+1))
        elif legend == 'potential':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_im_1.append("Z'' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_cycleno.append(str(np.round(np.average(self.df[i].E_avg), 2))+' V')


        if plot == 'w_data':
            fig = figure(figsize=(6, 8), dpi=120, facecolor='w', edgecolor='k')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(311, aspect='equal')
            ax1 = fig.add_subplot(312)
            ax2 = fig.add_subplot(313)
    
            colors = sns.color_palette("colorblind", n_colors=len(self.df))
            colors_real = sns.color_palette("Blues", n_colors=len(self.df)+2)
            colors_imag = sns.color_palette("Oranges", n_colors=len(self.df)+2)
    
            ### Nyquist Plot
            for i in range(len(self.df)):
                ax.plot(self.df[i].re, self.df[i].im, marker='o', ms=4, lw=2, color=colors[i], ls='-', alpha=.7, label=self.label_cycleno[i])
    
            ### Bode Plot
            if bode == 'on':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i+1], marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_re_1[i])
                    ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i+1], marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_im_1[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("Z', -Z'' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)

            elif bode == 're':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i+1], marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("Z' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)

            elif bode == 'log_re':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i+1], marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(Z') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)

            elif bode == 'im':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i+1], marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("-Z'' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)

            elif bode == 'log_im':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i+1], marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(-Z'') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)      

            elif bode == 'log':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i+1], marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_re_1[i])
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i+1], marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_im_1[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best', fontsize=10, frameon=False)

            ### Kramers-Kronig Relative Residuals    
            for i in range(len(self.df)):
                ax2.plot(np.log10(self.df[i].f), self.KK_rr_re[i]*100, color=colors_real[i+1], marker='D', ls='--', ms=6, alpha=.7, label=self.label_re_1[i])
                ax2.plot(np.log10(self.df[i].f), self.KK_rr_im[i]*100, color=colors_imag[i+1], marker='s', ls='--', ms=6, alpha=.7, label=self.label_im_1[i])
                ax2.set_xlabel("log(f) [Hz]")
                ax2.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")
                if legend == 'on' or legend == 'potential': 
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
            ax2.axhline(0, ls='--', c='k', alpha=.5)
            
            ### Setting ylims and write 'KK-Test' on RR subplot
            self.KK_rr_im_min = []
            self.KK_rr_im_max = []
            self.KK_rr_re_min = []
            self.KK_rr_re_max = []
            for i in range(len(self.df)):
                self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))    
            if np.min(self.KK_rr_im_min) > np.min(self.KK_rr_re_min):
                ax2.set_ylim(np.min(self.KK_rr_re_min)*100*1.5, np.max(np.abs(self.KK_rr_re_min))*100*1.5)
                ax2.annotate('Lin-KK', xy=[np.min(np.log10(self.df[0].f)), np.max(self.KK_rr_re_max)*100*.9], color='k', fontweight='bold')
            elif np.min(self.KK_rr_im_min) < np.min(self.KK_rr_re_min):
                ax2.set_ylim(np.min(self.KK_rr_im_min)*100*1.5, np.max(self.KK_rr_im_max)*100*1.5)
                ax2.annotate('Lin-KK', xy=[np.min(np.log10(self.df[0].f)), np.max(self.KK_rr_im_max)*100*.9], color='k', fontweight='bold')
                
            ### Figure specifics
            if legend == 'on' or legend == 'potential':
                ax.legend(loc='best', fontsize=10, frameon=False)
            ax.set_xlabel("Z' [$\Omega$]")
            ax.set_ylabel("-Z'' [$\Omega$]")
            if nyq_xlim != 'none':
                ax.set_xlim(nyq_xlim[0], nyq_xlim[1])
            if nyq_ylim != 'none':
                ax.set_ylim(nyq_ylim[0], nyq_ylim[1])
            #Save Figure
            if savefig != 'none':
                fig.savefig(savefig)

        ### Illustrating residuals only
    
        elif plot == 'residuals':
            colors = sns.color_palette("colorblind", n_colors=9)
            colors_real = sns.color_palette("Blues", n_colors=9)
            colors_imag = sns.color_palette("Oranges", n_colors=9)

            ### 1 Cycle
            if len(self.df) == 1:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax = fig.add_subplot(231)
                ax.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax.set_xlabel("log(f) [Hz]")
                ax.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")
                if legend == 'on' or legend == 'potential':
                    ax.legend(loc='best', fontsize=10, frameon=False)        
                ax.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and write 'KK-Test' on RR subplot
                self.KK_rr_im_min = np.min(self.KK_rr_im)
                self.KK_rr_im_max = np.max(self.KK_rr_im)
                self.KK_rr_re_min = np.min(self.KK_rr_re)
                self.KK_rr_re_max = np.max(self.KK_rr_re)
                if self.KK_rr_re_max > self.KK_rr_im_max:
                    self.KK_ymax = self.KK_rr_re_max
                else:
                    self.KK_ymax = self.KK_rr_im_max
                if self.KK_rr_re_min < self.KK_rr_im_min:
                    self.KK_ymin = self.KK_rr_re_min
                else:
                    self.KK_ymin = self.KK_rr_im_min
                if np.abs(self.KK_ymin) > self.KK_ymax:
                    ax.set_ylim(self.KK_ymin*100*1.5, np.abs(self.KK_ymin)*100*1.5)
                    if legend == 'on':
                        ax.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin)*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin)*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin) < self.KK_ymax:
                    ax.set_ylim(np.negative(self.KK_ymax)*100*1.5, np.abs(self.KK_ymax)*100*1.5)
                    if legend == 'on':
                        ax.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax*100*1.3], color='k', fontweight='bold')

                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 2 Cycles
            elif len(self.df) == 2:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                #cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])

                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 3 Cycles
            elif len(self.df) == 3:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.3], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 4 Cycles
            elif len(self.df) == 4:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(221)
                ax2 = fig.add_subplot(222)
                ax3 = fig.add_subplot(223)
                ax4 = fig.add_subplot(224)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                ax3.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best', fontsize=10, frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')

                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 5 Cycles
            elif len(self.df) == 5:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                ax4 = fig.add_subplot(234)
                ax5 = fig.add_subplot(235)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=18)
                ax4.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best', fontsize=10, frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best', fontsize=10, frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 6 Cycles
            elif len(self.df) == 6:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                ax4 = fig.add_subplot(234)
                ax5 = fig.add_subplot(235)
                ax6 = fig.add_subplot(236)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_xlabel("log(f) [Hz]")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best', fontsize=10, frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best', fontsize=10, frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best', fontsize=10, frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)
                          
            ### 7 Cycles
            elif len(self.df) == 7:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best', fontsize=10, frameon=False)
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best', fontsize=10, frameon=False)
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best', fontsize=10, frameon=False)
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_xlabel("log(f) [Hz]")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax7.legend(loc='best', fontsize=10, frameon=False)
                ax7.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)      
                           
            ### 8 Cycles
            elif len(self.df) == 8:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                ax8 = fig.add_subplot(338)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best', fontsize=10, frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best', fontsize=10, frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best', fontsize=10, frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_xlabel("log(f) [Hz]")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)                
                if legend == 'on' or legend == 'potential':
                    ax7.legend(loc='best', fontsize=10, frameon=False)        
                ax7.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 8
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_re[7]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_im[7]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax8.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax8.legend(loc='best', fontsize=10, frameon=False)        
                ax8.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[7]) > self.KK_ymax[7]:
                    ax8.set_ylim(self.KK_ymin[7]*100*1.5, np.abs(self.KK_ymin[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[7]) < self.KK_ymax[7]:
                    ax8.set_ylim(np.negative(self.KK_ymax[7])*100*1.5, np.abs(self.KK_ymax[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymax[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK, ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), self.KK_ymax[7]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 9 Cycles
            elif len(self.df) == 9:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                ax8 = fig.add_subplot(338)
                ax9 = fig.add_subplot(339)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on': 
                    ax1.legend(loc='best', fontsize=10, frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax2.legend(loc='best', fontsize=10, frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax3.legend(loc='best', fontsize=10, frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on': 
                    ax4.legend(loc='best', fontsize=10, frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax5.legend(loc='best', fontsize=10, frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax6.legend(loc='best', fontsize=10, frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                ax7.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax7.legend(loc='best', fontsize=10, frameon=False)
                ax7.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 8
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_re[7]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_im[7]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax8.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax8.legend(loc='best', fontsize=10, frameon=False)
                ax8.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 9
                ax9.plot(np.log10(self.df[8].f), self.KK_rr_re[8]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax9.plot(np.log10(self.df[8].f), self.KK_rr_im[8]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax9.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax9.legend(loc='best', fontsize=10, frameon=False)
                ax9.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[7]) > self.KK_ymax[7]:
                    ax8.set_ylim(self.KK_ymin[7]*100*1.5, np.abs(self.KK_ymin[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[7]) < self.KK_ymax[7]:
                    ax8.set_ylim(np.negative(self.KK_ymax[7])*100*1.5, np.abs(self.KK_ymax[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymax[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK, ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), self.KK_ymax[7]*100*1.2], color='k', fontweight='bold')

                if np.abs(self.KK_ymin[8]) > self.KK_ymax[8]:
                    ax9.set_ylim(self.KK_ymin[8]*100*1.5, np.abs(self.KK_ymin[8])*100*1.5)
                    if legend == 'on': 
                        ax9.annotate('Lin-KK, #9', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymin[8])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax9.annotate('Lin-KK ('+str(np.round(np.average(self.df[8].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymin[8])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[8]) < self.KK_ymax[8]:
                    ax9.set_ylim(np.negative(self.KK_ymax[8])*100*1.5, np.abs(self.KK_ymax[8])*100*1.5)
                    if legend == 'on': 
                        ax9.annotate('Lin-KK, #9', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymax[8])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax9.annotate('Lin-KK, ('+str(np.round(np.average(self.df[8].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[8].f)), self.KK_ymax[8]*100*1.2], color='k', fontweight='bold')  
                        
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)
            else:
                print('Too many spectras, cannot plot all. Maximum spectras allowed = 9')




def leastsq_errorfunc(params, w, re, im, circuit, weight_func):
    '''
    Sum of squares error function for the complex non-linear least-squares fitting procedure (CNLS). The fitting function (lmfit) will use this function to iterate over
    until the total sum of errors is minimized.
    
    During the minimization the fit is weighed, and currently three different weigh options are avaliable:
        - modulus
        - unity
        - proportional
    
    Modulus is generially recommended as random errors and a bias can exist in the experimental data.
        
    Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)

    Inputs
    ------------
    - params: parameters needed for CNLS
    - re: real impedance
    - im: Imaginary impedance
    - circuit:
      The avaliable circuits are shown below, and this this parameter needs it as a string.
        - C
        - Q
        - R-C
        - R-Q
        - RC
        - RQ
        - R-RQ
        - R-RQ-RQ
        - R-RQ-Q
        - R-(Q(RW))
        - R-(Q(RM))
        - R-RC-C
        - R-RC-Q
        - R-RQ-Q
        - R-RQ-C
        - RC-RC-ZD
        - R-TLsQ
        - R-RQ-TLsQ
        - R-TLs
        - R-RQ-TLs
        - R-TLQ
        - R-RQ-TLQ
        - R-TL
        - R-RQ-TL
        - R-TL1Dsolid (reactive interface with 1D solid-state diffusion)
        - R-RQ-TL1Dsolid

    - weight_func
      Weight function
        - modulus
        - unity
        - proportional
    '''
    if circuit == 'C':
        re_fit = elem_C_fit(params, w).real
        im_fit = -elem_C_fit(params, w).imag
    elif circuit == 'Q':
        re_fit = elem_Q_fit(params, w).real
        im_fit = -elem_Q_fit(params, w).imag
    elif circuit == 'R-C':
        re_fit = cir_RsC_fit(params, w).real
        im_fit = -cir_RsC_fit(params, w).imag
    elif circuit == 'R-Q':
        re_fit = cir_RsQ_fit(params, w).real
        im_fit = -cir_RsQ_fit(params, w).imag
    elif circuit == 'RC':
        re_fit = cir_RC_fit(params, w).real
        im_fit = -cir_RC_fit(params, w).imag
    elif circuit == 'RQ':
        re_fit = cir_RQ_fit(params, w).real
        im_fit = -cir_RQ_fit(params, w).imag
    elif circuit == 'R-RQ':
        re_fit = cir_RsRQ_fit(params, w).real
        im_fit = -cir_RsRQ_fit(params, w).imag
    elif circuit == 'R-RQ-RQ':
        re_fit = cir_RsRQRQ_fit(params, w).real
        im_fit = -cir_RsRQRQ_fit(params, w).imag
    elif circuit == 'R-RC-C':
        re_fit = cir_RsRCC_fit(params, w).real
        im_fit = -cir_RsRCC_fit(params, w).imag
    elif circuit == 'R-RC-Q':
        re_fit = cir_RsRCQ_fit(params, w).real
        im_fit = -cir_RsRCQ_fit(params, w).imag
    elif circuit == 'R-RQ-Q':
        re_fit = cir_RsRQQ_fit(params, w).real
        im_fit = -cir_RsRQQ_fit(params, w).imag
    elif circuit == 'R-RQ-C':
        re_fit = cir_RsRQC_fit(params, w).real
        im_fit = -cir_RsRQC_fit(params, w).imag
    elif circuit == 'R-(Q(RW))':
        re_fit = cir_Randles_simplified_Fit(params, w).real
        im_fit = -cir_Randles_simplified_Fit(params, w).imag
    elif circuit == 'R-(Q(RM))':
        re_fit = cir_Randles_uelectrode_fit(params, w).real
        im_fit = -cir_Randles_uelectrode_fit(params, w).imag
    elif circuit == 'C-RC-C':
        re_fit = cir_C_RC_C_fit(params, w).real
        im_fit = -cir_C_RC_C_fit(params, w).imag
    elif circuit == 'Q-RQ-Q':
        re_fit = cir_Q_RQ_Q_Fit(params, w).real
        im_fit = -cir_Q_RQ_Q_Fit(params, w).imag
    elif circuit == 'RC-RC-ZD':
        re_fit = cir_RCRCZD_fit(params, w).real
        im_fit = -cir_RCRCZD_fit(params, w).imag
    elif circuit == 'R-TLsQ':
        re_fit = cir_RsTLsQ_fit(params, w).real
        im_fit = -cir_RsTLsQ_fit(params, w).imag
    elif circuit == 'R-RQ-TLsQ':
        re_fit = cir_RsRQTLsQ_Fit(params, w).real
        im_fit = -cir_RsRQTLsQ_Fit(params, w).imag
    elif circuit == 'R-TLs':
        re_fit = cir_RsTLs_Fit(params, w).real
        im_fit = -cir_RsTLs_Fit(params, w).imag
    elif circuit == 'R-RQ-TLs':
        re_fit = cir_RsRQTLs_Fit(params, w).real
        im_fit = -cir_RsRQTLs_Fit(params, w).imag
    elif circuit == 'R-TLQ':
        re_fit = cir_RsTLQ_fit(params, w).real
        im_fit = -cir_RsTLQ_fit(params, w).imag
    elif circuit == 'R-RQ-TLQ':
        re_fit = cir_RsRQTLQ_fit(params, w).real
        im_fit = -cir_RsRQTLQ_fit(params, w).imag
    elif circuit == 'R-TL':
        re_fit = cir_RsTL_Fit(params, w).real
        im_fit = -cir_RsTL_Fit(params, w).imag
    elif circuit == 'R-RQ-TL':
        re_fit = cir_RsRQTL_fit(params, w).real
        im_fit = -cir_RsRQTL_fit(params, w).imag
    elif circuit == 'R-TL1Dsolid':
        re_fit = cir_RsTL_1Dsolid_fit(params, w).real
        im_fit = -cir_RsTL_1Dsolid_fit(params, w).imag
    elif circuit == 'R-RQ-TL1Dsolid':
        re_fit = cir_RsRQTL_1Dsolid_fit(params, w).real
        im_fit = -cir_RsRQTL_1Dsolid_fit(params, w).imag
    else:
        print('Circuit is not defined in leastsq_errorfunc()')
        
    error = [(re-re_fit)**2, (im-im_fit)**2] #sum of squares
    
    #Different Weighing options, see Lasia
    if weight_func == 'modulus':
        weight = [1/((re_fit**2 + im_fit**2)**(1/2)), 1/((re_fit**2 + im_fit**2)**(1/2))]
    elif weight_func == 'proportional':
        weight = [1/(re_fit**2), 1/(im_fit**2)]
    elif weight_func == 'unity':
        unity_1s = []
        for k in range(len(re)):
            unity_1s.append(1) #makes an array of [1]'s, so that the weighing is == 1 * sum of squres.
        weight = [unity_1s, unity_1s]
    else:
        print('weight not defined in leastsq_errorfunc()')
        
    S = np.array(weight) * error #weighted sum of squares 
    return S




#IMPORT THE DATA FILE IN THE FORM OF AN MPT FILE
#working on adjusting to mpt if not an mpt file to begin with
def importer(path, data, mask_front, mask_back):
    mpt = mpt_data(path, data, mask = [mask_front, mask_back])
    df = mpt.df_raw
    mpt.mpt_plot()
    return [mpt, df]


def mpt_plot(mpt):
    mpt.EIS_plot()


def guess(guess_package):
    
    #SINGLE ITERATION OF THE GUESS PROCESS
    #USE THIS FUNCTION TO GET CLOSER TO THE IDEAL COEFFICIENTS FOR Rs, R, n, fs, R2, n2, fs2
    #REPEAT THIS FUNCTION UNTIL THE THRESHOLD IS ACHEIVED
    
    params = Parameters()
    
    #adding to the parameters package to send to the fitting function
    params.add('Rs', value=guess_package[0], min=guess_package[0]*.01, max=guess_package[0]*100)
    params.add('R', value=guess_package[1], min=guess_package[1]*.1, max=guess_package[1]*10)
    params.add('n', value=guess_package[2], min=.65, max=1.2)
    params.add('fs', value=guess_package[3], min=10**0.5, max=10**6)
    params.add('R2', value=guess_package[4], min=guess_package[4]*.1, max=guess_package[4]*10)
    params.add('n2', value=guess_package[5], min=.65, max=1.2)
    params.add('fs2', value=guess_package[6], min=10**-2, max=10**1)
    
    #Call to the fitting function given by PyEIS
    mpt_data.EIS_fit(params=params, circuit='R-RQ-RQ', weight_func='modulus')
    
    #maybe take a look at the plots,may help for accuracy, don't really need it...
    #mpt_data.EIS_plot(fitting = 'on')
    
    
    #print out the values
    print(mpt_data.fit_Rs)
    print()
    print(mpt_data.fit_R)
    print(mpt_data.fit_n)
    print(mpt_data.fit_fs)
    print()
    print(mpt_data.fit_R2)
    print(mpt_data.fit_n2)
    print(mpt_data.fit_fs2)
    
    #export the new guess package
    guess_package =  ([mpt_data.fit_Rs[0],mpt_data.fit_R[0],mpt_data.fit_n[0],mpt_data.fit_fs[0],mpt_data.fit_R2[0],mpt_data.fit_n2[0],mpt_data.fit_fs2[0]])
    return guess_package



#THIS VERIFIES WHETHER OR NOT WE'VE ACHEIVED A SATISFACTORY COEFFICIENT PACKAGE
#IF THIS DOESN'T RETURN TRUE, WE RUN THE GUESSER UNTIL IT DOES
def thresh_verif(before, after):
    try:
        total = 0
        for i in range(len(before)):
            total += (before[i] - after[i])
        print(total)    
        return abs(total) <= 1e-10
    except IndexError as e:
        #IF LISTS AREN'T THE SAME LENGTH
        print("Lists are not the same length")
        return



#ITERATIVE GUESSER
def guesser(Rs_guess,R_guess,n_guess,fs_guess,R2_guess,n2_guess,fs2_guess):
    guess_package = [Rs_guess, R_guess, n_guess, fs_guess, R2_guess, n2_guess, fs2_guess]
    new_guess = guess(guess_package)
    while not thresh_verif(guess_package, new_guess):
        guess_package = new_guess
        new_guess = guess(new_guess)
        print(new_guess)
    return new_guess