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
import sys, traceback
pd.options.mode.chained_assignment = None
import statistics as stat
from os import listdir
from os.path import isfile, join


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


pd.options.mode.chained_assignment = None

#TAKEN FROM PYEIS LIBRARY
def extract_mpt(path, EIS_name):
    '''
    Extracting PEIS and GEIS data files from EC-lab '.mpt' format, coloums are renames following correct_text_EIS()
    
    Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)
    '''
    EIS_init = pd.read_csv(path+EIS_name, sep='\t', nrows=1,header=0,names=['err'], encoding='latin1') #findes line that states skiplines
    EIS_test_header_names = pd.read_csv(path+EIS_name, sep='\t', skiprows=int(EIS_init.err[0][18:-1])-1, encoding='latin1') #locates number of skiplines
    names_EIS = []
    for j in range(len(EIS_test_header_names.columns)):
        names_EIS.append(correct_text_EIS(EIS_test_header_names.columns[j])) #reads coloumn text
    return pd.read_csv(path+EIS_name, sep='\t', skiprows=int(EIS_init.err[0][18:-1]), names=names_EIS, encoding='latin1')

#TAKEN FROM PYEIS LIBRARY
def correct_text_EIS(text_header):
    '''Corrects the text of '*.mpt' and '*.dta' files into readable parameters without spaces, ., or /
    
    <E_we> = averaged Wew value for each frequency
    <I> = Averaged I values for each frequency
    |E_we| = module of Ewe
    |I_we| = module of Iwe
    Cs/F = Capacitance caluculated using an R+C (series) equivalent circuit
    Cp/F = Capacitance caluculated using an R-C (parallel) equivalent circuit
    Ref.:
        - EC-Lab User's Manual
    
    Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)
    '''
    if text_header == 'freq/Hz' or text_header == '  Freq(Hz)':
        return 'f'
    elif text_header == 'Re(Z)/Ohm' or text_header == "Z'(a)":
        return 're'
    elif text_header == '-Im(Z)/Ohm' or text_header == "Z''(b)":
        return 'im'
    elif text_header == '|Z|/Ohm':
        return 'Z_mag'
    elif text_header == 'Phase(Z)/deg':
        return 'Z_phase'
    elif text_header == 'time/s' or text_header == 'Time(Sec)':
        return 'times'
    elif text_header == '<Ewe>/V' or text_header == 'Bias':
        return 'E_avg'
    elif text_header == '<I>/mA':
        return 'I_avg'
    elif text_header == 'Cs/F':
        return 'Cs' ####
    elif text_header == 'Cp/F':
        return 'Cp'
    elif text_header == 'cycle number':
        return 'cycle_number'
    elif text_header == 'Re(Y)/Ohm-1':
        return 'Y_re'
    elif text_header == 'Im(Y)/Ohm-1':
        return 'Y_im'
    elif text_header == '|Y|/Ohm-1':
        return 'Y_mag'
    elif text_header == 'Phase(Y)/deg':
        return 'Y_phase'
    elif text_header == 'Time':
        return 'times'
    elif text_header == 'Freq':
        return 'f'
    elif text_header == 'Zreal':
        return 're'
    elif text_header == 'Zimag':
        return 'im'
    elif text_header == 'Zmod':
        return 'Z_mag'
    elif text_header == 'Vdc':
        return 'E_avg'
    elif text_header == 'Idc':
        return 'I_avg'
    elif text_header == 'I/mA':
        return 'ImA'
    elif text_header == 'Ewe/V':
        return 'EweV'
    elif text_header == 'half cycle':
        return 'half_cycle'
    elif text_header == 'Ns changes':
        return 'Ns_changes'
    else:
        return text_header


class mpt_data:
    def __init__(self, path, data, cycle='off', mask=['none','none'], gph_width = 6.4, gph_height = 4.8):
        self.path = path
        self.data = data
        self.width = gph_width
        self.height = gph_height
        self.df_raw0 = []
        self.cycleno = []
        self.mask = mask
        self.counter = 0
        self.low_error = 0
        for j in range(len(data)):
            if data[j].find(".mpt") != -1: #file is a .mpt file
                self.df_raw0.append(extract_mpt(path=path, EIS_name=data[j])) #reads all datafiles
            else:
                print('Data file(s) could not be identified')

            self.cycleno.append(self.df_raw0[j].cycle_number)
            if np.min(self.cycleno[j]) <= np.max(self.cycleno[j-1]):
                if j > 0: #corrects cycle_number except for the first data file
                    self.df_raw0[j].update({'cycle_number': self.cycleno[j]+np.max(self.cycleno[j-1])}) #corrects cycle number

        self.df_raw = [i for i in self.df_raw0][0]
        self.df_raw = self.df_raw.assign(w = 2*np.pi*self.df_raw.f)

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
    
    #DEFINE SIZE OF GRAPH
    #BE CAREFUL AS THIS DOESN'T DETERMINE THE DIMENSIONS OF THE WINDOW
    #THIS ONLY SETS THE DIMENSIONS OF THE SIZE OF THE WINDOW
    #ACTUAL GRAPH DATA DIMENSIONS CAN BE ADJUSTED IN THE PLOTTING FUNCTION
    def set_gph_width(self, new_width):
        self.width = new_width
        return
    
    def set_gph_height(self, new_height):
        self.height = new_height
        return

    def set_new_gph_dims(self, new_width, new_height):
        self.set_gph_width(new_width)
        self.set_gph_height(new_height)
        return
    
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

    def window_masker(self, x_window, y_window):
        adj_re = self.df_raw[(self.df_raw['re']<x_window[1]) & (self.df_raw['re']>x_window[0])]
        adj_mpt = adj_re[(adj_re['im']<y_window[1]) & (adj_re['im']>y_window[0])]
        return [max(adj_mpt['f']), min(adj_mpt['f'])]
    
    #PLOTTING FUNCTION
    def mpt_plot(self, fitting='off', rr='off', legend='on', x_window = 'none', y_window = 'none'):
        
        #Figure Initialization
        
        fig = plt.figure(dpi=120, figsize = [self.width, self.height], facecolor='w', edgecolor='w')
        fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
        ax = fig.add_subplot(211, aspect='equal')
        
        
        ### Figure specifics
        if legend == 'on': 
            ax.legend(loc='best', fontsize=10, frameon=False)
        ax.set_xlabel("Z' [$\Omega$]")
        ax.set_ylabel("-Z'' [$\Omega$]")
        if x_window != 'none':
            ax.set_xlim(x_window[0], x_window[1])
        if y_window != 'none':
            ax.set_ylim(y_window[0], y_window[1])
        
        #Color initialization
        colors = sns.color_palette("colorblind", n_colors=len(self.df))
        colors_real = sns.color_palette("Blues", n_colors=len(self.df)+2)
        colors_imag = sns.color_palette("Oranges", n_colors=len(self.df)+2)
    
        #Label functions
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

        ### Relative Residuals on Fit
        """
        if rr=='on':
            ax2 = fig.add_subplot(212)
            if fitting == 'off':
                print('Fitting has not been performed, thus the relative residuals cannot be determined')
            elif fitting == 'on':
                self.rr_real = []
                self.rr_imag = []
                for i in range(len(self.df)):
                    self.rr_real.append(residual_real(re=self.df[i].re.values, fit_re=self.circuit_fit[i].real, fit_im=-self.circuit_fit[i].imag))
                    self.rr_imag.append(residual_imag(im=self.df[i].im.values, fit_re=self.circuit_fit[i].real, fit_im=-self.circuit_fit[i].imag))
                    if legend == 'on':
                        ax2.plot(np.log10(self.df[i].f), self.rr_real[i]*100, color=colors_real[i], marker='D', ms=6, lw=1, ls='--', label='#'+str(i+1))
                        ax2.plot(np.log10(self.df[i].f), self.rr_imag[i]*100, color=colors_imag[i], marker='s', ms=6, lw=1, ls='--',label='')
                    elif legend == 'potential':
                        ax2.plot(np.log10(self.df[i].f), self.rr_real[i]*100, color=colors_real[i], marker='D', ms=6, lw=1, ls='--', label=str(np.round(np.average(self.df[i].E_avg.values),2))+' V')
                        ax2.plot(np.log10(self.df[i].f), self.rr_imag[i]*100, color=colors_imag[i], marker='s', ms=6, lw=1, ls='--',label='')

                    ax2.axhline(0, ls='--', c='k', alpha=.5)
                    ax2.set_xlabel("log(f) [Hz]")
                    ax2.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")

                #Automatic y-limits limits
                self.rr_im_min = []
                self.rr_im_max = []
                self.rr_re_min = []
                for i in range(len(self.df)): # needs to be within a loop if cycles have different number of data points     
                    self.rr_im_min = np.min(self.rr_imag[i])
                    self.rr_im_max = np.max(self.rr_imag[i])
                    self.rr_re_min = np.min(self.rr_real[i])
                    self.rr_re_max = np.max(self.rr_real[i])
                if self.rr_re_max > self.rr_im_max:
                    self.rr_ymax = self.rr_re_max
                else:
                    self.rr_ymax = self.rr_im_max
                if self.rr_re_min < self.rr_im_min:
                    self.rr_ymin = self.rr_re_min
                else:
                    self.rr_ymin  = self.rr_im_min
                if np.abs(self.rr_ymin) > np.abs(self.rr_ymax):
                    ax2.set_ylim(self.rr_ymin *100*1.5, np.abs(self.rr_ymin)*100*1.5)
                    ax2.annotate("$\Delta$Z'", xy=(np.log10(np.min(self.df[0].f)), np.abs(self.rr_ymin )*100*1.2), color=colors_real[-1], fontsize=12)
                    ax2.annotate("$\Delta$-Z''", xy=(np.log10(np.min(self.df[0].f)), np.abs(self.rr_ymin )*100*0.9), color=colors_imag[-1], fontsize=12)
                elif np.abs(self.rr_ymin) < np.abs(self.rr_ymax):
                    ax2.set_ylim(np.negative(self.rr_ymax)*100*1.5, np.abs(self.rr_ymax)*100*1.5)                    
                    ax2.annotate("$\Delta$Z'", xy=(np.log10(np.min(self.df[0].f)), np.abs(self.rr_ymax)*100*1.2), color=colors_real[-1], fontsize=12)
                    ax2.annotate("$\Delta$-Z''", xy=(np.log10(np.min(self.df[0].f)), np.abs(self.rr_ymax)*100*0.9), color=colors_imag[-1], fontsize=12)
    
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best', fontsize=10, frameon=False)
        """


        ### Nyquist Plot
        ax.plot(self.df[0].re, self.df[0].im, marker='o', ms=4, lw=2, color=colors[i], ls='-', label=self.label_cycleno[i])
        if fitting == 'on':
            real = []
            imag = []
            for i in self.circuit_fit[0]:
                #print(i.real)
                real.append(i.real)
                #print(i.imag)
                imag.append(-i.imag)
            ax.plot(real, imag, lw=0, marker='o', ms=8, mec='r', mew=1, mfc='none', label='')
        plt.show()


    #FITTING THE FREQUENCY ONTO THE GRAPH. FLIP SWITCH ON PLOT FUNCT TO DISPLAY
    def mpt_fit(self, params, circuit, weight_func='modulus', nan_policy='raise'):
        assert circuit == 'R-RQ-RQ'

        self.Fit = []
        self.circuit_fit = []
        self.fit_E = []
        self.fit_Rs = []
        self.fit_R = []
        self.fit_n = []
        self.fit_R2 = []
        self.fit_n2 = []
        self.fit_fs = []
        self.fit_fs2 = []
        self.fit_Q = []
        self.fit_Q2 = []
        self.fit_Q3 = []
        self.fit_n3 = []
        self.fit_fs3 = []

        for i in range(len(self.df)):
            self.Fit.append(minimize(self.leastsq_errorfunc, params, method='leastsq', args=(self.df[i].w.values, self.df[i].re.values, self.df[i].im.values, circuit, weight_func), nan_policy=nan_policy, maxfev=9999990))
            print(report_fit(self.Fit[i]))
            #print(self.Fit)
            #self.fit_E.append(np.average(self.df[i].E_avg))
        
        
        for i in range(len(self.df)):
            if "'fs'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()) and "'fs3'" in str(self.Fit[i].params.keys()):
                #print('HERE')
                self.circuit_fit.append(cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value, R2=self.Fit[i].params.get('R2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value, Q3='none', n3=self.Fit[i].params.get('n3').value,fs3 = self.Fit[i].params.get('fs3').value)),
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R.append(self.Fit[i].params.get('R').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
                self.fit_fs.append(self.Fit[i].params.get('fs').value)
                self.fit_R2.append(self.Fit[i].params.get('R2').value)
                self.fit_n2.append(self.Fit[i].params.get('n2').value)
                self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                self.fit_fs3.append(self.Fit[i].params.get('fs3').value)
                self.fit_Q.append(1/(self.fit_R[0] * (self.fit_fs[0] * 2 * np.pi)**self.fit_n[0])) 
                self.fit_Q2.append(1/(self.fit_R2[0] * (self.fit_fs2[0] * 2 * np.pi)**self.fit_n2[0])) 
                self.fit_n3.append(self.Fit[i].params.get('n3').value) 
                #print(self.fit_fs3[0] * 2 * np.pi)
                #print(self.fit_n3[0])
                self.fit_Q3.append(1/((self.fit_fs3[0] * 2 * np.pi)**self.fit_n3[0])) 
            else:
                print("Circuit Error, check inputs")
                break
        #print(self.circuit_fit)

    def leastsq_errorfunc(self, params, w, re, im, circuit, weight_func):
        re_fit = cir_RsRQRQ_fit(params, w).real
        im_fit = -cir_RsRQRQ_fit(params, w).imag
        error = [(re-re_fit)**2, (im-im_fit)**2] #sum of squares
        print('MPT FILE : ', self.data[0], ' ERROR: ', sum(error))
        self.low_error = sum(error)
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
    
    #Updated Guesser
    def guesser(self, csv_container = None):
        Rs_guess = 1e3
        R_guess = 1 
        n_guess = 0.8 
        fs_guess = 1 
        R2_guess = 100 
        n2_guess = 0.8 
        fs2_guess = 0.2 
        n3_guess = 0.8
        fs3_guess = 1


        params = Parameters()
        #adding to the parameters package to send to the fitting function
        params.add('Rs', value=Rs_guess, min=Rs_guess*.01, max=10**6)
        params.add('R', value=R_guess, min=Rs_guess*.1, max=10**6)
        params.add('n', value=n_guess, min=.65, max=1)
        params.add('fs', value=fs_guess, min=10**0.5, max=10**6)
        params.add('R2', value=R2_guess, min=R2_guess*.1, max=10**6)
        params.add('n2', value=n2_guess, min=.65, max=1)
        params.add('fs2', value=fs2_guess, min=10**-2, max=10**6)
        params.add('n3', value=n3_guess, min=.65, max=1)
        params.add('fs3', value=fs3_guess, min=10**-2, max=10**6)
        self.mpt_fit(params, circuit = 'R-RQ-RQ')

        counter = 0

        while self.low_error >= 10000 and counter <= 100:        
            try:
                counter += 1
                print('ITERATION NO. : ', counter)
                Rs_guess = self.fit_Rs[0]

                R_guess = self.fit_R[0]
                n_guess = self.fit_n[0]
                fs_guess = self.fit_fs[0]

                R2_guess = self.fit_R2[0]
                n2_guess = self.fit_n2[0]
                fs2_guess = self.fit_fs2[0]

                n3_guess = self.fit_n3[0]
                fs3_guess = self.fit_fs3[0]

                guess_package = [Rs_guess, R_guess, n_guess, fs_guess, R2_guess, n2_guess, fs2_guess, n3_guess, fs3_guess]
                #adding to the parameters package to send to the fitting function
                params = Parameters()
                params.add('Rs', value=guess_package[0], min=guess_package[0]*.01, max=guess_package[0]*100)
                params.add('R', value=guess_package[1], min=guess_package[1]*.1, max=guess_package[1]*10)
                params.add('n', value=guess_package[2], min=.65, max=1)
                params.add('fs', value=guess_package[3], min=10**0.5, max=10**6)
                params.add('R2', value=guess_package[4], min=guess_package[4]*.1, max=guess_package[4]*10)
                params.add('n2', value=guess_package[5], min=.65, max=1)
                params.add('fs2', value=guess_package[6], min=10**-2, max=10**1)
                params.add('n3', value=guess_package[7], min=.65, max=1)
                params.add('fs3', value=guess_package[8], min=10**-2, max=10**1)
                self.mpt_fit(params, circuit = 'R-RQ-RQ')


            except KeyboardInterrupt:
                print('Interrupted!!')
                #print([self.fit_Rs[0],self.fit_R[0],self.fit_n[0],self.fit_Q[0],self.fit_R2[0],self.fit_n2[0],self.fit_Q2[0]])
        #self.set_new_gph_dims(50,50)
        #self.mpt_plot(fitting = 'on')
        self.fitted = pd.DataFrame({'file':self.data,
                    'fit_R':self.fit_Rs,
                "fit_Rs":self.fit_R,
                "fit_n":self.fit_n,
                "fit_Q":self.fit_Q,
                "fit_R2":self.fit_R2,
                "fit_n2":self.fit_n2,
                "fit_Q2":self.fit_Q2,
                "fit_n3":self.fit_n3,
                "fit_Q3":self.fit_Q3})
        out_name = 'fitted_' + self.data[0][:-4]
        if csv_container:
            self.fitted.to_csv(csv_container+out_name, sep='\t')
            return self.fitted
        return self.fitted


def cir_RsRQRQ_fit(params, w):
    '''
    Fit Function: -Rs-RQ-RQ-
    Return the impedance of an Rs-RQ circuit. See details under cir_RsRQRQ()
    
    Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
    '''
    if str(params.keys())[10:].find("'R'") == -1: #if R == 'none':
        Q = params['Q']
        n = params['n']
        fs = params['fs']
        R = (1/(Q*(2*np.pi*fs)**n))
    if str(params.keys())[10:].find("'Q'") == -1: #elif Q == 'none':
        R = params['R']
        n = params['n']
        fs = params['fs']
        Q = (1/(R*(2*np.pi*fs)**n))
    if str(params.keys())[10:].find("'n'") == -1: #elif n == 'none':
        R = params['R']
        Q = params['Q']
        fs = params['fs']
        n = np.log(Q*R)/np.log(1/(2*np.pi*fs))
    if str(params.keys())[10:].find("'fs'") == -1: #elif fs == 'none':
        R = params['R']
        Q = params['Q']
        n = params['n']

    if str(params.keys())[10:].find("'R2'") == -1: #if R == 'none':
        Q2 = params['Q2']
        n2 = params['n2']
        fs2 = params['fs2']
        R2 = (1/(Q2*(2*np.pi*fs2)**n2))
    if str(params.keys())[10:].find("'Q2'") == -1: #elif Q == 'none':
        R2 = params['R2']
        n2 = params['n2']
        fs2 = params['fs2']
        Q2 = (1/(R2*(2*np.pi*fs2)**n2))
    if str(params.keys())[10:].find("'n2'") == -1: #elif n == 'none':
        R2 = params['R2']
        Q2 = params['Q2']
        fs2 = params['fs2']
        n2 = np.log(Q2*R2)/np.log(1/(2*np.pi*fs2))
    if str(params.keys())[10:].find("'fs2'") == -1: #elif fs == 'none':
        R2 = params['R2']
        Q2 = params['Q2']
        n2 = params['n2']
    
    n3 = params['n3']
    Q3 = (1/((2*np.pi*params['fs3'])**n3))
    Rs = params['Rs']
    return Rs + (R/(1+R*Q*(w*1j)**n)) + (R2/(1+R2*Q2*(w*1j)**n2)) + (1/(Q3*(w*1j))**n3)


def cir_RsRQRQ(w, Rs, R='none', Q='none', n='none', fs='none', R2='none', Q2='none', n2='none', fs2='none', Q3 = 'none', fs3 = 'none', n3 = 'none'):

    
    if Q == 'none':
        Q = (1/(R*(2*np.pi*fs)**n))
    
   
    if Q2 == 'none':
        Q2 = (1/(R2*(2*np.pi*fs2)**n2))
    

    if Q3 == 'none':
        Q3 = (1/(1*(2*np.pi*fs3)**n3))
    
        
    return Rs + (R/(1+R*Q*(w*1j)**n)) + (R2/(1+R2*Q2*(w*1j)**n2)) + (1/(Q3*(w*1j))**n3)


def full_graphing(path, lst = None):
    bad_mpts = []
    if not lst:
        path_files = [f for f in listdir(path) if isfile(join(path, f)) if f[-3:] == 'mpt']
        for i in path_files:
            try:
                print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                masked_mpt.set_new_gph_dims(30,30)
                masked_mpt.mpt_plot()
                plt.show()
            except ValueError:
                bad_mpts.append(i)
                print(i, ' was a bad file, could not find a mask')
        if bad_mpts:
            print(bad_mpts, " are a list of bad mpts. You may want to take a closer look at them")
    if type(lst) == list:
        for i in lst:
            try:
                print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                masked_mpt.set_new_gph_dims(30,30)
                masked_mpt.mpt_plot()
            except ValueError:
                bad_mpts.append(i)
                print(i, ' was a bad file, could not find a mask')
            if bad_mpts:
                print(bad_mpts, " are a list of bad mpts. You may want to take a closer look at them")
                

def the_ringer(path, single_file):
        print('WHOLE THING')
        ex_mpt = mpt_data(path,[single_file])
        ex_mpt.mpt_plot()
        
        print('FAST MASK')
        print(ex_mpt.fast_mask())
        fast_masked_mpt = mpt_data(path, [single_file], mask = ex_mpt.fast_mask())
        fast_masked_mpt.mpt_plot()
        
        print('MASKER0')
        print(ex_mpt.masker0())
        masker0_mpt = mpt_data(path, [single_file], mask = ex_mpt.masker0())
        masker0_mpt.mpt_plot()
        
        print('MASKER')
        print(ex_mpt.masker())
        masker_mpt = mpt_data(path, [single_file], mask = ex_mpt.masker())
        masker_mpt.mpt_plot()
                       
#PATH takes in a string that leads to the files
#CSV_CONTAINER takes an additional path that leads to a separate folder which will contain all the fitted coefficients
#if you want to just fit a single mpt or a list of mpts, you can use LST for specific fittings
#TAKE_CSV for when you want to export a csv

def auto_fit(path, entry, csv_container = None):
    bad_mpts = []

    fitteds = []
    if type(entry) == list:
        for i in entry:
            try:
                #print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                fitteds.append(masked_mpt.guesser(csv_container = csv_container))
            except ValueError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                bad_mpts.append(i)
                ex_mpt.mpt_plot()
                print(i, ' was a bad file, could not find a mask')
            except TypeError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                ex_mpt.guesser(csv_container = csv_container)
                print(i, ' was fittable, but could not obtain a mask')
    
    if type(entry) == str:
        path_files = [f for f in listdir(path) if isfile(join(path, f)) if f[-3:] == 'mpt']
        for i in path_files:
            try:
                #print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                fitteds.append(masked_mpt.guesser(csv_container = csv_container))
            except ValueError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                bad_mpts.append(i)
                ex_mpt.mpt_plot()
                print(i, ' was a bad file, could not find a mask')
            except TypeError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                ex_mpt.guesser(csv_container = csv_container)
                print(i, ' was fittable, but could not obtain a mask')
    
    to_export = pd.concat(fitteds)
    return to_export

def path_listing(path):
    path_files = [f for f in listdir(path) if isfile(join(path, f)) if f[-3:] == 'mpt']
    for i in path_files:
        print(i)