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



#IMPORT THE DATA FILE IN THE FORM OF AN MPT FILE
#working on adjusting to mpt if not an mpt file to begin with
def importer(path, data, mask_front, mask_back):
    mpt_data = EIS_exp(path, data, mask = [mask_front, mask_back])
    df = mpt_data.df_raw
    mpt_data.EIS_plot()
    return [mpt_data, df]


def mpt_plot(mpt):
    mpt_data.EIS_plot()