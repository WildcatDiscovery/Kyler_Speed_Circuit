#ASSISTING FUNCTIONS FOR THE IMPEDANCE_DATA IPYTHON NOTEBOOK
#DERIVED FROM KRISTIAN KNUDSEN'S PYEIS REPO
#HELPING TO FIT POINTS FROM A NYQUIST PLOT IN THE FORM OF A MPT FILE

#IMPORT NESSECARY LIBRARIES
from PyEIS import *


class EIS_exp:
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
        if len(self.df_raw0) == 1:
            self.df_raw = self.df_raw0[0]
        elif len(self.df_raw0) == 2:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1]], axis=0)
        elif len(self.df_raw0) == 3:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2]], axis=0)
        elif len(self.df_raw0) == 4:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3]], axis=0)
        elif len(self.df_raw0) == 5:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4]], axis=0)
        elif len(self.df_raw0) == 6:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5]], axis=0)
        elif len(self.df_raw0) == 7:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6]], axis=0)
        elif len(self.df_raw0) == 8:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7]], axis=0)
        elif len(self.df_raw0) == 9:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8]], axis=0)
        elif len(self.df_raw0) == 10:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9]], axis=0)
        elif len(self.df_raw0) == 11:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10]], axis=0)
        elif len(self.df_raw0) == 12:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], axis=0)
        elif len(self.df_raw0) == 13:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11], self.df_raw0[12]], axis=0)
        elif len(self.df_raw0) == 14:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], self.df_raw0[12], self.df_raw0[13], axis=0)
        elif len(self.df_raw0) == 15:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], self.df_raw0[12], self.df_raw0[13], self.df_raw0[14], axis=0)
        else:
            print("Too many data files || 15 allowed")
        self.df_raw = self.df_raw.assign(w = 2*np.pi*self.df_raw.f) #creats a new coloumn with the angular frequency

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

    def EIS_plot(self, bode='off', fitting='off', rr='off', nyq_xlim='none', nyq_ylim='none', legend='on', savefig='none'):
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

    def EIS_fit(self, params, circuit, weight_func='modulus', nan_policy='raise'):
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



#IMPORT THE DATA FILE IN THE FORM OF AN MPT FILE
#working on adjusting to mpt if not an mpt file to begin with
def importer(path, data, mask_front, mask_back):
    mpt_data = EIS_exp(path, data, mask = [mask_front, mask_back])
    df = mpt_data.df_raw
    mpt_data.EIS_plot()
    return [mpt_data, df]


def mpt_plot(mpt):
    mpt.EIS_plot()

