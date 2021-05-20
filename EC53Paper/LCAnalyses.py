def AutoCorrelation(dates, pfluxes, noises,interval=30,normalize=1,revolving=1,Full=0,Resample=0,Custom=0):
    import numpy as np
    
    
    x = np.linspace(0,dates[-1],dates[-1])
    resampling_function = np.arange(dates[0],dates[-1],interval)
    if Custom:
        resampling_function = Resample
#    print(resampling_function)
    dates = np.array(dates)
    pfluxes = np.array(pfluxes)
    noises = np.array(noises)
    resampled_fluxes = []
    resampled_fluxes.append(pfluxes[0])
    resampled_noises = []
    resampled_noises.append(noises[0])

    for point in resampling_function[1:]:
        bridging_points_x = [dates[np.where(dates < point)][-1],dates[np.where(dates > point)][0]]
        bridging_points_y = [pfluxes[np.where(dates < point)][-1],pfluxes[np.where(dates > point)][0]]
        a = [bridging_points_x[0],bridging_points_y[0]]
        b = [bridging_points_x[1],bridging_points_y[1]]
        slope = (b[1]-a[1])/(b[0]-a[0])
        resampled_flux = a[1] + slope*(point-a[0])
        noise_cal = (noises[np.where(dates < point)][-1]**2+noises[np.where(dates > point)][0]**2)
        resampled_noise = np.sqrt(noises[np.where(dates < point)][-1]**2 + (slope**2)*noise_cal)
        resampled_fluxes.append(resampled_flux)
        resampled_noises.append(resampled_noise)

    
    r_dates = np.array(resampling_function)
    r_fluxes = np.array(resampled_fluxes)
    r_noises = np.array(resampled_noises)
    
    ACF = []
    mean_r_fluxes = np.mean(r_fluxes)

    if normalize:
        r_fluxes = [x-mean_r_fluxes for x in r_fluxes]
    ACF.append(np.correlate(r_fluxes,r_fluxes))
    if normalize:
        ACF[0] = ACF[0]/np.sum([x**2 for x in r_fluxes])
    
    for j in range(1,int(len(r_fluxes)/2)):
        i = len(r_fluxes)-j
        
        shifted = list(r_fluxes)[i:]+list(r_fluxes)[:i]
        point_ACF = np.correlate(r_fluxes,np.array(shifted))
        if not revolving:
            point_ACF = np.correlate(r_fluxes[j:],r_fluxes[:-j])
        if normalize:
            point_ACF = point_ACF/np.sum([x**2 for x in r_fluxes])
        ACF.append(point_ACF)
    
    ACF = np.array(ACF)
    
    return(resampling_function,resampled_fluxes,resampled_noises, ACF)

def LCResample(dates, pfluxes, noises,r_func):
    import numpy as np
    
    #x = np.linspace(0,dates[-1],dates[-1])
    resampling_function = r_func
    dates = np.array(dates)
    pfluxes = np.array(pfluxes)
    noises = np.array(noises)
    resampled_fluxes = []
    resampled_noises = []

    for point in resampling_function[0:]:
        bridging_points_x = [dates[np.where(dates < point)][-1],dates[np.where(dates > point)][0]]
        bridging_points_y = [pfluxes[np.where(dates < point)][-1],pfluxes[np.where(dates > point)][0]]
        a = [bridging_points_x[0],bridging_points_y[0]]
        b = [bridging_points_x[1],bridging_points_y[1]]
        slope = (b[1]-a[1])/(b[0]-a[0])
        resampled_flux = a[1] + slope*(point-a[0])
        noise_cal = (noises[np.where(dates < point)][-1]**2+noises[np.where(dates > point)][0]**2)
        resampled_noise = np.sqrt(noises[np.where(dates < point)][-1]**2 + (slope**2)*noise_cal)
        resampled_fluxes.append(resampled_flux)
        resampled_noises.append(resampled_noise)
    
    r_dates = np.array(resampling_function)
    r_fluxes = np.array(resampled_fluxes)
    r_noises = np.array(resampled_noises)

    return(r_dates, r_fluxes, r_noises)

def LCSEDflux(dates,fluxes,errors,t_dates,wl):
    import numpy as np
    import LCAnalyses as LCA
    import scipy.constants as cons
    
    wlsplit = wl.split("_")
    if wlsplit[1] == "J":
        wvs = 1.25
    if wlsplit[1] == "H":
        wvs = 1.65
    if wlsplit[1] == "K":
        wvs = 2.2
    if wlsplit[1] == "W1":
        wvs = 3.4
    if wlsplit[1] == "W2":
        wvs = 4.6
    
#    t_dates, r_fluxes, r_noises, ACF = LCA.AutoCorrelation(dates,fluxes,errors,Custom=1, Resample=t_dates)
    t_dates, r_fluxes, r_noises = LCA.LCResample(dates,fluxes,errors,t_dates)
    SED_fluxes = [x/wvs*(cons.c*1e+6) for x in r_fluxes]
    SED_errors = [x/wvs*(cons.c*1e+6) for x in r_noises]
    
    return(t_dates, SED_fluxes,SED_errors)



def LCStringlength(dates, fluxes, noises, periods,wl = '850', Phaseplot=0,progress=0):
    import numpy as np
    import matplotlib.pyplot as plt

    delta_dates = []
    string_lengths = []
    string_length_vars = []

    for i in range(0,len(periods)):
        coupled = [[dates[0],fluxes[0],noises[0]]]
        phased_dates = [x%periods[i] for x in dates]
    
        for j in range(0,len(phased_dates)):
            coupled.append([phased_dates[j],fluxes[j],noises[j]])
        coupled.sort(key=lambda coupled: coupled[0])
        s_fluxes = [x[1] for x in coupled]
        s_phase = [x[0]/periods[i] for x in coupled]
        s_noises = [x[2] for x in coupled]
        string_length = 0
        string_length_var = 0

        for k in range(1,len(s_fluxes)):
            delta_phase = s_phase[k] - s_phase[k-1]
            delta_flux = s_fluxes[k] - s_fluxes[k-1]
            delta_flux_error = np.sqrt(s_noises[k]**2 + s_noises[k-1]**2)
            line_segment = np.sqrt(delta_phase**2 + delta_flux**2)
            string_length = string_length + line_segment
            string_length_var = string_length_var + delta_flux_error**2
            if progress:
                if k == int(1*len(s_fluxes)/5):
                    print("tying up 20%")
                if k == int(2*len(s_fluxes)/5):
                    print("tying up 40%")
                if k == int(3*len(s_fluxes)/5):
                    print("tying up 60%")
                if k == int(4*len(s_fluxes)/5):
                    print("tying up 80%")
        string_lengths.append(string_length)
        string_length_vars.append(string_length_var)
        if periods[i] == Phaseplot:
            plt.errorbar(s_phase, s_fluxes,yerr=s_noises,fmt='o-',color='green',ecolor='green')
            #plt.plot(s_phase,s_fluxes)
            plt.xlabel("Phase")
            plt.ylabel("Flux [Jy/beam]")
            if wl != '850':
                plt.ylabel("Mangitude [mag]")
                plt.gca().invert_yaxis()
        #print("DONE!")
    return(string_lengths,np.sqrt(string_length_vars))

def LCDiffCoeff(dates, pfluxes, noises):
    import numpy as np
    
    dates = np.array(dates)
    pfluxes = np.array(pfluxes)
    noises = np.array(noises)
    d_dates = []
    d_fluxes = []
    d_noises = []

    for i in range(1,len(dates)):
        d_dates.append((dates[i]-dates[i-1])/2+dates[i-1])
        d_fluxes.append((pfluxes[i]-pfluxes[i-1])/(dates[i]-dates[i-1]))
        d_noises.append(np.sqrt(noises[i]**2+noises[i-1]**2)/(dates[i]-dates[i-1]))
        
    return(d_dates,d_fluxes,d_noises)

from astropy.timeseries import LombScargle as LS
import numpy as np

class LCLombscargle():
    
    def __init__(self,dates,pfluxes,noises,peak=0,window=0):

        self.dates = dates
        self.pfluxes = pfluxes
        self.noises = noises
        self.ls_pf = LS(dates, pfluxes, noises)
        freqs,powers = self.ls_pf.autopower()
        self.freqs = np.array(freqs)
        self.powers = np.array(powers)
        self.periods = 1.0/self.freqs
        self.residue_times = 0

        self.detected_freq = freqs[np.where(powers == powers.max())][0]
        self.FAPs = self.ls_pf.false_alarm_probability(self.powers)
        self.FAP_single = self.ls_pf.false_alarm_probability(np.max(self.powers),method='single')
        self.FAP1 = self.ls_pf.false_alarm_level(0.01)
        self.FAP001 = self.ls_pf.false_alarm_level(0.00001)
        if peak != 0:
            #print(self.periods)
            power_cut1 = np.where(self.periods>peak[0])[0][-1]
            power_cut2 = np.where(self.periods<peak[1])[0][0]
            #print(power_cut1, power_cut2)
            power_cut = self.powers[power_cut2:power_cut1+1]
            RPeak = power_cut.max()
            detected_freq_new = self.freqs[np.where(self.powers == RPeak)][0]
            self.detected_freq = detected_freq_new
            print("Current Set Peak : " + str(1/self.detected_freq) +" days")
        self.FAP_peak = self.ls_pf.false_alarm_probability(self.powers[np.where(self.freqs == self.detected_freq)])

    
        
    def Periodogram(self):
        #if nterms == 1:
        #    FAP1 = ls_pf.false_alarm_level(0.01)
        #    FAP01 = ls_pf.false_alarm_level(0.001)
        #    FAPs = ls_pf.false_alarm_probability(power)
        return(self.freqs, self.powers)
    
    def mulsin(self,x, *sinpar_in_ori):
        # input must be [offset, freq1, amp1, iniphase1 (, freq2, amp2, iniphase2, .....)]
        y = np.array([0.0 for t in x])
        if len(sinpar_in_ori) == 1:
            sinpar_in = sinpar_in_ori[0]
        else:
            sinpar_in = sinpar_in_ori

        if np.size(sinpar_in)%4 == 0:
            sinpar = [0]
            for i in range(0,len(sinpar_in),4):
                sinpar[0] += sinpar_in[i]
                sinpar.append(sinpar_in[i+1])
                sinpar.append(sinpar_in[i+2])
                sinpar.append(sinpar_in[i+3])
        else:
            sinpar = sinpar_in
        
        y +=  sinpar[0]
        #print(sinpar)
        for i in range(1,len(sinpar),3):
            y += np.array([sinpar[i+1]*np.sin(2*np.pi*sinpar[i]*(x)-sinpar[i+2]) for x in x])
        #print(len(y))
        return(y)
        

    def modelpara(self):
        
        model_raw_para = self.ls_pf.model_parameters(self.detected_freq)
        #print(model_raw_para)
        amp = np.sqrt(model_raw_para[1]**2 + model_raw_para[2]**2)
        ini_phase = np.arctan(model_raw_para[2]/model_raw_para[1])
        #print(ini_phase)
        if model_raw_para[2] > 0 and model_raw_para[1] < 0:
            ini_phase += np.pi
            #print("case 1")
            
        if model_raw_para[2] < 0 and model_raw_para[1] < 0:
            ini_phase -= np.pi
            #print("case 2")
            
        ini_phase= -ini_phase
        output_para = [model_raw_para[0]+np.mean(self.pfluxes),self.detected_freq,amp,ini_phase]
        sinpar = output_para
        #output : [offset,Freqeuncy, amp, ini_phase]
        return(sinpar)
    
#    def residue(self):
#        currentmodel =self.ls_pf.model(self.dates,self.detected_freq)
#        residue = pfluxes - currentmodel
#        self.residue_times += 1
#        return(dates, residue)
        
    def window(self):
        ls_pf = LS(self.dates, [1 for x in self.pfluxes], self.noises,fit_mean=False, center_data=False)
        freqs, powers = ls_pf.autopower()
        return(freqs,powers)
        
    def sinfit(self,inipar_in,bounds=0):
        import scipy.optimize as op
        #output = []

        offset = 0
        inipar = [0]
        #print("here")
        #print(inipar_in[0])
        for i in range(0,len(inipar_in),4):
            inipar[0] += inipar_in[i]
            inipar.append(inipar_in[i+1])
            inipar.append(inipar_in[i+2])
            inipar.append(inipar_in[i+3])
        offset = inipar[0]
            # reorder the input parameters to [offset, P1, A1, Pi1, P2, A2, Pi2.....]
        print(inipar)
        if bounds:
            sol_para, sol_cov = op.curve_fit(self.mulsin, self.dates, self.pfluxes, sigma=self.noises, p0=inipar,maxfev=3000, bounds=bounds)
            sol_err = np.sqrt(np.diag(sol_cov))
            return(sol_para,sol_err)
        sol_para, sol_cov = op.curve_fit(self.mulsin, self.dates, self.pfluxes, sigma=self.noises, p0=inipar,absolute_sigma=True,maxfev=3000)
        sol_err = np.sqrt(np.diag(sol_cov))
        #print(len(sol_para))
        #for i in range(0,len(sol_para),4):
        #    output.append(sol_para[i:i+4])
        return(sol_para, sol_err)
    

    
def LCSampling(r_dates, r_fluxes, r_noises, rate=0.8):
    import random as rd
    import numpy as np
    indices = [0]
    a = 0
    while a < rate:
        idx = int(rd.uniform(1,len(r_dates)))
        if np.sum([x == idx for x in indices]):
            continue
        indices.append(idx)
        a = len(indices)/len(r_dates)

    dates = np.array(r_dates)[indices]
    fluxes = np.array(r_fluxes)[indices]
    noises = np.array(r_noises)[indices]
    coupled = []
    for j in range(0,len(dates)):
        coupled.append([dates[j], fluxes[j], noises[j]])
    coupled.sort(key = lambda coupled: coupled[0])
    dates = [x[0] for x in coupled]
    fluxes = [x[1] for x in coupled]
    noises = [x[2] for x in coupled]
    
    
    return(dates, fluxes, noises)
    
def LCLinearfit(dates, fluxes, noises, xnoises=[0],plot=0,mod_date='',region='',k='',nonoise=0):
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.optimize as op
    import statistics as stat
    import LCPlot as LCP
    from scipy import odr
    
    def func(x,a,b):
        return(a*x + b)
    def odrfunc(p,x):
        a,b = p
        return(a*x + b)
    #sol_para = [1,1]
    #sol_err
    if not xnoises[0]:
        sol_para, sol_cov = op.curve_fit(func,dates,fluxes,sigma=noises,absolute_sigma=True)
    else:
        func_model = odr.Model(odrfunc)
        data = odr.RealData(dates, fluxes ,sx=xnoises, sy=noises)
        odrc = odr.ODR(data, func_model, beta0=[0.,1.])
        out = odrc.run()
        
        return(out.beta, np.diag(out.cov_beta))
    
    if nonoise:
        sol_para, sol_cov = op.curve_fit(func,dates,fluxes)
        
    sol_err = np.sqrt(np.diag(sol_cov))
    
    if plot:
        x = np.linspace(min(dates)-100,max(dates) + 100)
        LCP.LCplot(dates[:], fluxes[:], noises[:], mod_date, region, k, wl='850', mag=0, option='', save=0, newplot=1, color='green',ecolor = 'black',stringcolor='green',markersize=5,capsize=3,label="",
          legend=0,title=0)
        plt.plot(x,func(x,sol_para[0],sol_para[1]))
        plt.xlim(min(dates)-100,max(dates)+100)
        plt.plot()
        plt.xticks(rotation=20)
        plt.title('Light curve ')
        #plt.annotate('S$_{fit}$'+' ='+' %4.2f'%(sol_para[0]),xy=[0.15,0.82], xycoords="figure fraction")
        plt.annotate('S$_{fit}$'+'/'+'$\Delta$S$_{fit}$ ='+' %4.2f'%(sol_para[0]/sol_err[0]),xy=[0.15,0.77], xycoords="figure fraction")
    
    return(sol_para, sol_err)

def Gaussfit(x, y, yn=0):
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.optimize as op
    import statistics as stat
    import LCPlot as LCP
    from scipy import asarray as ar,exp
    n = len(x)
    mean = np.mean(x)
    sigma = np.std(x)
    
    def gaus(x,a,x0,sigma):
        return a*exp(-(x-x0)**2/(2*sigma**2))
    
    if yn:
        popt, pcov = op.curve_fit(gaus, x, y, p0=[np.max(y),mean,sigma], sigma=yn, absolute_sigma=True)
    else:
        popt, pcov = op.curve_fit(gaus, x, y, p0=[np.max(y),mean,sigma])
    
    return(popt, np.sqrt(np.diag(pcov)))
    
    
def LCSamplingFunction(region, mod_date,dates=[],fluxes=[],noises=[], noise=0.02,wl='850', option='', hist=1,given=0):
    import numpy as np
    import LCAnalyses as LCA
    import importlib
    import numpy as np
    import LCAnalyses as LCA
    import LCCall as LCC
    import LCPlot as LCP
    import matplotlib.pyplot as plt
    import random
    importlib.reload(LCA)
    importlib.reload(LCC)
    importlib.reload(LCP)

    datadir = "/mnt/Secdrive/TSmain/"+region
    output_dir = datadir + '/Results_'+mod_date+'/'
    
    if not given:
        datadir = "/mnt/Secdrive/TSmain/"+region
        output_dir = datadir + '/Results_'+mod_date+'/'
        snames, dic = LCC.JCMTTransient(region,mod_date)
        k = snames[0]
        dates, fluxes, noises = LCC.JCMTTRansientSource(k,region,wl,mod_date)
        dates =[x-dates[0] for x in dates]
        fluxes = [1.0 for x in fluxes]
        noises = [noise for x in noises]
        if wl == '850':
            noises = [0.02 for x in noises]
        if wl == '450':
            noises = [0.05 for x in noises]
    if given:
        noises = [noise for x in noises]
        fluxes = [1.0 for x in fluxes]
    detected_periods = []
    for i in range(0,1000):
        fluxes = [x + random.gauss(0,y) for x, y in zip(fluxes,noises)]
        freqs, powers, faps, fap1, fap10 = LCA.LCLombscargle(dates, fluxes, noises)
        periods = 1/freqs
        detected_periods.append(float(periods[np.where(powers == powers.max())]))
        if i == 0:
            plt.close()
            plt.plot(periods,powers)
            plt.plot([min(periods)-1,max(periods)+1],[fap1,fap1],color='black',)
            plt.xlim([min(periods),max(periods)])
            plt.xscale('log')
            plt.savefig(output_dir + region +'_obsf_PS_'+wl+'.pdf')
            plt.close()
        
    plt.close()
    plt.hist(detected_periods,alpha=0.3,color='green',bins=30,label='Peak position of Periodogram')
    plt.legend()

    #plt.yscale('log')
    plt.xlabel('Period')
    plt.ylabel('N(peak)')
    #plt.xscale('log')
    plt.xticks()
    #plt.annotate()

    plt.savefig(output_dir+region+'_ObsF_PShist_'+mod_date+'_'+wl+'.pdf')
    print(output_dir+region+'_ObsF_PShist_'+mod_date+'_'+wl+'.pdf')

    return(0)

def LCPeriodogramerr(freqs, power, peakpos = 0,plot=0):
    import scipy.optimize as op
    import numpy as np
    import matplotlib.pyplot as plt
    
    #peakpos0 = 1/peakpos
    peakpos1 = 1/(peakpos-100)
    peakpos2 = 1/(peakpos+100)
    print(peakpos1, peakpos2)
    def gaussian(x, pos, amp, sig,base = 0):
        return(amp*np.exp(-np.power(x - pos, 2.) / (2*np.power(sig, 2.)))+base)
    power_cut1 = np.where(freqs>peakpos1)[0][0]
    power_cut2 = np.where(freqs<peakpos2)[0][-1]
    print(power_cut1,power_cut2)
    power_cut = power[power_cut2:power_cut1+1]
    RPeak = np.max(power_cut)
    #RPeak = power[np.where(periods == peakpos)]
    RPeakpos = np.where(power == RPeak)[0][0]
    pos_Aroundpeak = freqs[(RPeakpos-2):(RPeakpos+3)]
    Aroundpeak = power[(RPeakpos-2):(RPeakpos+3)]
    initial_guess = [pos_Aroundpeak[2],RPeak,0.001]
    sol_para, sol_cov = op.curve_fit(gaussian,pos_Aroundpeak,Aroundpeak,p0=initial_guess)
    sol_err = np.sqrt(np.diag(sol_cov))
    err_range = np.array([0,0,0])
    err_range[0] =np.round(1/(sol_para[2]+pos_Aroundpeak[2]))
    err_range[2]=np.round(1/(pos_Aroundpeak[2]-sol_para[2]))
    err_range[1]=np.round(1/(pos_Aroundpeak[2]))
    print(RPeakpos,err_range-err_range[1])
    return(sol_para,sol_err,err_range)
    
def LCTimeWinCut(dates, fluxes, noises, initial_date,final_date):
    import numpy as np
    
    index = np.where(dates > initial_date)[0][0]
    index2= np.where(dates < final_date)[0][-1]
    print(dates[index],dates[index2])
    dates_cut = [x-initial_date for x in dates[index:(index2+1)]]
    
    fluxes_cut = fluxes[index:(index2+1)]
    noises_cut = noises[index:(index2+1)]
    
    return(np.array(dates_cut), np.array(fluxes_cut), np.array(noises_cut))
    
import numpy as np
class LCphaseup:
    def __init__(self,dates,data,wl, offset,offset2=0,offset_date=150,period=530,cutter = 2457640,cutter2=2458760):
        self.offset = offset
        self.offset2 = offset2
        self.wl = wl
        self.dates = np.array(dates)
        self.data = data
        
        if wl.split("_")[0] == 'UKIRT':
            #self.cut_ind=  55  # to 56 for first period
            self.cut_ind = np.where(self.dates<cutter)[0][-1]+1
        if wl.split("_")[0] == "SCUBA2":
            #self.cut_ind= 7 # 7 for first period
            self.cut_ind = np.where(self.dates<cutter)[0][-1]+1
            self.data = -2.5*np.log10(self.data/self.data[0])
            if wl.split("_")[1] == "850":
                self.data = self.data*5.5
            if wl.split("_")[1] == "450":
                self.cut_ind = np.where(self.dates<cutter)[0][-1]+1
                self.data = self.data*3
        if wl.split("_")[0] == "WISE":
            self.data = self.data[:-3]
            self.dates = self.dates[:-3]
            #self.cut_ind = 13
            self.cut_ind = np.where(self.dates<cutter)[0][-1]+1
            
        self.cut_ind2 = np.where(self.dates<cutter2)[0][-1]+1
        self.xplot = (self.dates+offset_date)%period/period
        self.yplot = list(self.data[:self.cut_ind]-(self.offset))+list(self.data[self.cut_ind:])
        #self.cut_ind2 = len(self.dates)-1
        if offset2:
            self.yplot = list(self.data[:self.cut_ind]-(self.offset)) + list(self.data[self.cut_ind:self.cut_ind2]) +list(self.data[self.cut_ind2:]-(self.offset2))
        self.ylabel = self.wl.split("_")[1] + ' [mag]'
        
    def CutLC(self):
        xout = self.xplot[:self.cut_ind]
        yout = self.yplot[:self.cut_ind]
        return(xout,yout)
    
    def Phasepickup(self,phases):
        xout= np.array(self.xplot)[np.where((self.xplot>phases[0]) & (self.xplot<phases[1]))]
        yout= np.array(self.yplot)[np.where((self.xplot>phases[0]) & (self.xplot<phases[1]))]
        return(xout,yout)

    
    def Phasesplit(self,phases,show = 0):
        output = []
        outputx = []
        for i in range(1,len(phases)):
            output.append(np.array(self.yplot)[np.where((self.xplot>=phases[i-1]) & (self.xplot<phases[i]))])
            outputx.append(np.array(self.xplot)[np.where((self.xplot>=phases[i-1]) & (self.xplot<phases[i]))])
            if show:
                print(len(output[-1]))
        return(outputx,output)

    def Phasextend(self):
        extendx = list(self.xplot-1) + list(self.xplot) + list(self.xplot+1)
        self.xplott = extendx
        extendy = self.yplot +self.yplot + self.yplot
        self.yplott = extendy
        return(len(self.xplott),len(self.yplott))


