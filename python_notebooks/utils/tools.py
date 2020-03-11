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


from utils.data_extraction import *
from utils.lin_kk import *

#print(path, data)
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
        re_mid, im_mid  = mean(skeleton['re']), mean(skeleton['im'])
        a = skeleton[abs(skeleton['re']) <= re_mid * .5]
        b = skeleton[abs(skeleton['im']) <= im_mid * .5]
        c = pd.concat([a, b]).drop_duplicates()
        return [c['f'].max(), c['f'].min()]

    def kk_masker(self,number = 1):

        num_RC='auto' 
        legend='on'
        plot='residuals'
        bode='off'
        nyq_xlim='none'
        nyq_ylim='none'
        weight_func='Boukamp'



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

        self.KK_circuit_fit = []
        self.KK_rr_re = []
        self.KK_rr_im = []
        functs = []
        for i in range(2,81):
            functs.append('KK_RC'+str(i))

        for i in range(len(self.df)):
            cir_num = int(self.number_RC[i])
            cir_funct = eval(functs[cir_num - 2])
            self.KK_circuit_fit.append(cir_funct(w=self.df[0].w, Rs=self.Lin_KK_Fit[0].params.get('Rs').value, R_values=self.KK_R[0], t_values=self.t_const[0]))
            if cir_num >= 81:
                print('RC simulation circuit not defined')
                print('   Number of RC = ', self.number_RC)
            self.KK_rr_re.append(residual_real(re=self.df[i].re, fit_re=self.KK_circuit_fit[i].real, fit_im=-self.KK_circuit_fit[i].imag)) #relative residuals for the real part
            self.KK_rr_im.append(residual_imag(im=self.df[i].im, fit_re=self.KK_circuit_fit[i].real, fit_im=-self.KK_circuit_fit[i].imag)) #relative residuals for the imag part



        self.kk_df = pd.DataFrame({'f':np.log10(self.df_raw.f), 're':self.KK_rr_re[0]*100, 'im':self.KK_rr_im[0]*100})
        self.kk_df['difference'] = abs(self.kk_df['re'] - self.kk_df['im'])
        diff_mean = self.kk_df['difference'].mean()
        masked_df = self.kk_df[self.kk_df['difference'] < diff_mean * number]
        return [10**masked_df['f'].max(),10**masked_df['f'].min()]

    def masker0(self):
        skeleton = self.df_raw.iloc[:,0:3]
        re_lim, im_lim  = max(skeleton['re']) * .6, max(skeleton['im'] * .6)
        a = skeleton[(skeleton['re']) <= re_lim]
        b = skeleton[(skeleton['im']) <= im_lim]
        c = pd.concat([a, b]).drop_duplicates()

        return [max(c['f']), min(c['f'])]

    def masker(self, num_bins = 5):

        c = self.df_raw.iloc[:,0:3]
        for cols in c.columns.tolist()[1:]:
            c = c.ix[c[cols] > 0]

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

    
    #PLOTTING FUNCTION
    def mpt_plot(self, fitting='off', rr='off', legend='on', x_window = 'none', y_window = 'none'):
        
        #Figure Initialization
        fig = figure(dpi=120, figsize = [self.width, self.height], facecolor='w', edgecolor='w')
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



        ### Nyquist Plot
        for i in range(len(self.df)):
            ax.plot(self.df[i].re, self.df[i].im, marker='o', ms=4, lw=2, color=colors[i], ls='-', label=self.label_cycleno[i])
            if fitting == 'on':
                ax.plot(self.circuit_fit[i].real, -self.circuit_fit[i].imag, lw=0, marker='o', ms=8, mec='r', mew=1, mfc='none', label='')
        
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
    def guesser(self, csv_container, to_csv = False):
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

        while self.low_error >= 100000 and counter <= 100:        
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
        self.set_new_gph_dims(50,50)
        self.mpt_plot(fitting = 'on')
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
        if to_csv == True:
            self.fitted.to_csv(csv_container+out_name, sep='\t')
            return self.fitted
        return self.fitted


"""    #DETERMINE THE OPTIMAL MASK THROUGH LINEAR KRAMER KRONIG ANALYSIS      
    def Lin_KK(self, num_RC='auto', legend='on', plot='residuals', bode='off', nyq_xlim='none', nyq_ylim='none', weight_func='Boukamp', savefig='none'):
        #NEED TO REDOCUMENT
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
"""

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
        
        print('KK MASKER')
        print(kk_masker(ex_mpt))
        kk_masked_mpt = mpt_data(path, [single_file], mask = kk_masker(ex_mpt))
        kk_masked_mpt.mpt_plot()
        
        print('FAST MASK')
        print(fast_mask(ex_mpt))
        fast_masked_mpt = mpt_data(path, [single_file], mask = fast_mask(ex_mpt))
        fast_masked_mpt.mpt_plot()
        
        print('MASKER0')
        print(masker0(ex_mpt))
        masker0_mpt = mpt_data(path, [single_file], mask = masker0(ex_mpt))
        masker0_mpt.mpt_plot()
        
        print('MASKER')
        print(masker(ex_mpt))
        masker_mpt = mpt_data(path, [single_file], mask = masker(ex_mpt))
        masker_mpt.mpt_plot()
                       
#PATH takes in a string that leads to the files
#CSV_CONTAINER takes an additional path that leads to a separate folder which will contain all the fitted coefficients
#if you want to just fit a single mpt or a list of mpts, you can use LST for specific fittings
#TAKE_CSV for when you want to export a csv

def auto_fit(path, csv_container, lst = None, take_csv = False):
    bad_mpts = []
    fitted_files = [f for f in listdir(csv_container) if isfile(join(csv_container, f)) if f[:9] == 'fitted_DE']
    path_files = [f for f in listdir(path) if isfile(join(path, f)) if f[-3:] == 'mpt']
    if not lst:
        for i in path_files:
            try:
                #print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name not in fitted_files:
                    masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                    masked_mpt.guesser(csv_container = csv_container, to_csv = take_csv)
                else:
                    print(i, ' has already been fitted!!')
                    continue
            except ValueError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name in fitted_files:
                    print(i, ' has already been fitted!!')
                    continue
                bad_mpts.append(i)
                ex_mpt.mpt_plot()
                print(i, ' was a bad file, could not find a mask')
            except TypeError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name in fitted_files:
                    print(i, ' has already been fitted!!')
                    continue
                ex_mpt.guesser(csv_container = csv_container, to_csv = take_csv)
                print(i, ' was fittable, but could not obtain a mask')
    for i in lst:
            try:
                #print(i, ' was a permissible file')
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name not in fitted_files:
                    masked_mpt = mpt_data(path,[i], mask = ex_mpt.masker())
                    masked_mpt.guesser(csv_container = csv_container, to_csv = take_csv)
                else:
                    print(i, ' has already been fitted!!')
                    continue
            except ValueError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name in fitted_files:
                    print(i, ' has already been fitted!!')
                    continue
                bad_mpts.append(i)
                ex_mpt.mpt_plot()
                print(i, ' was a bad file, could not find a mask')
            except TypeError:
                ex_mpt = mpt_data(path,[i])
                out_name = 'fitted_' + ex_mpt.data[0][:-4]
                if out_name in fitted_files:
                    print(i, ' has already been fitted!!')
                    continue
                ex_mpt.guesser(csv_container = csv_container, to_csv = take_csv)
                print(i, ' was fittable, but could not obtain a mask')