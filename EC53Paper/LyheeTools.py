

def Shuffle(deck, order=0):
    import random as rd
    import numpy as np
    deck_shuffled = np.asarray(deck)
    rd.shuffle(deck_shuffled)
    
    #if order:
    #    if len(order) != len(deck) :
    #        raise TypeError("Input order should have the same size with the deck")
    return(list(deck_shuffled))
    

def JD2date(mjd, fmt='mjd'):
    import julian
#    dates = []
    dt = julian.from_jd(mjd, fmt='mjd')
    if fmt == 'jd':
        dt = julian.from_jd(mjd, fmt='jd')
        
    date = str('{:04d}').format(dt.year) + str('{:02d}').format(dt.month) +str('{:02d}').format(dt.day)
    return(date)
    

def TCMultiLS(k,wl,dates,pfluxes,noises,output_dir,errsize=3,marksize=5,close_up=5,y_reverse=0):
    import importlib
    import numpy as np
    from astropy.stats import LombScargle as LS
    import matplotlib.pyplot as plt
    from datetime import date
    import statistics as stat
    import TCTriggerFunctions as TF

    obs_pfluxes = tuple(pfluxes)
    ls_pf = LS(dates, pfluxes, noises)
    freq, power = ls_pf.autopower()
    
    detected = 1
    comp_index = 0
    detection_threshold  = 0.01
    detection_level = ls_pf.false_alarm_level(detection_threshold)
    detection_FAP = ls_pf.false_alarm_probability(power.max())
    detected_periods = []
    fits = []
    grid = np.linspace(-100,2000,3000)
    print(detection_FAP)
    while power.max() > detection_level:
        
        detection_FAP = ls_pf.false_alarm_probability(power.max())

        detected = 1
        comp_index = comp_index + 1
        print(str(comp_index) + " component FAP: " + str(detection_FAP))
        detected_freq = float(freq[np.where(power == power.max())])
        detected_period = 1/detected_freq
        print(str(comp_index) + " component period: " +str(detected_period))
        fit = ls_pf.model(grid,detected_freq)
        fit_in_dates = ls_pf.model(dates,detected_freq)
        fits.append(fit)
        detected_periods.append(detected_period)
        
        plt.figure()
        plt.errorbar([x for x in dates],pfluxes,yerr=noises,fmt='o',color='g',ecolor='black',capsize=errsize, markersize= marksize)
        plt.plot(grid, fit, color = 'Darkturquoise')
        plt.annotate("~" + "{:5.0F}".format(detected_period)+ " days",xy=[0.15,0.77], xycoords="figure fraction") 
        plt.xlim([-10,np.array(dates).max()+10])
        plt.xlabel('Relative Dates [Days]')
        plt.ylabel('Peak Flux [Jy/beam]')
        plt.title('Light Curve of ' + k + '(' + wl +')'+" comp. #" + str(comp_index))
        if y_reverse == 1:
            plt.gca().invert_yaxis()
            plt.ylabel("Magnitudes")
        plt.savefig(output_dir  + k +'_'+ 'LC_'+ wl  +'_'+ str(comp_index) + 'comp.pdf')
        plt.close()

        plt.figure()
        p_tough_thr = ls_pf.false_alarm_level(0.001)
        p_thr = detection_level
        plt.plot(freq*365.24, power,color='g')
        plt.plot([-1,freq.max()*365.24],[p_thr, p_thr], color='black')
        plt.plot([-1,freq.max()*365.24],[p_tough_thr, p_tough_thr], color='gray')
        #if not(detected_freq == freq.min()):
        plt.text(freq[np.where(power == power.max())], power.max()+0.05, "~" + "{:5.0F}".format(detected_period)+ " days")
        plt.xlabel("Cycles per year")
        plt.ylabel("Power")
        plt.title("Power Spectrum of " + k + ' (' + wl +')'+ " comp. #" + str(comp_index),color='black')

        plt.xlim([0,freq.max()*365.24])
        plt.ylim([0,1])
        plt.savefig(output_dir+ k + '_PS_' + wl +'_'+str(comp_index) + 'comp.pdf')
        plt.xlim([detected_freq*365.24-close_up/2,detected_freq*365.24+close_up/2])
        if detected_freq < close_up:
            plt.xlim([0,close_up])
        plt.savefig(output_dir+ k + '_PS_' + wl +'_'+str(comp_index) + 'comp_peak.pdf')
        plt.close()

        pfluxes = np.sum([pfluxes,[-x for x in fit_in_dates]],axis=0)
        ls_pf = LS(dates, pfluxes, noises)
        freq, power = ls_pf.autopower()
        
        final_fit = [x for x in np.sum(fits,axis = 0)]
        plt.plot(grid,final_fit,color='Darkturquoise')
        plt.xlim([-10,np.array(dates).max()+10])
        plt.xlabel("Relative Dates [Days]")
        plt.ylabel("Peak Flux [Jy/beam]")
        if y_reverse == 1:
            plt.gca().invert_yaxis()
            plt.ylabel("Magnitudes")
        plt.savefig(output_dir + k +'_LC_'+ wl + '_allcomp.pdf')
        plt.close()
    
    plt.figure()    
    plt.errorbar([x for x in dates],obs_pfluxes,yerr=noises,fmt='o',color='g',ecolor='black',capsize=errsize, markersize= marksize)
    transparency = np.linspace(0,1,comp_index+1)
    #for i in range(0,comp_index):
    #    plt.plot(grid,fits[i][:],color="Darkturquoise",alpha = transparency[i])
    final_fit = [x for x in np.sum(fits,axis = 0)]
    plt.plot(grid,final_fit,color='Darkturquoise')
    plt.xlim([-10,np.array(dates).max()+10])
    plt.xlabel("Relative Dates [Days]")
    plt.ylabel("Peak Flux [Jy/beam]")
    if y_reverse == 1:
        plt.gca().invert_yaxis()
        plt.ylabel("Magnitudes")
    plt.savefig(output_dir + k +'_LC_'+ wl + '_allcomp.pdf')
    plt.close()
    
    print("Done!")
    print("Detected components :" + str(comp_index))
    return(grid, fits,comp_index)


def TCPhaseFold(source_name,region,pfluxes,noises,grid,fitted_sine,Ydates,detected_freq,datadir,filenameoption,highlight=0,y_reverse=0,errsize=3,symsize=5):
    import matplotlib.pyplot as plt
    import numpy as np
    
    detected_period = 1/detected_freq
    Ydates = np.array(Ydates)
    Obsphases = int(Ydates.max()*detected_freq)+1
    resultpwd = datadir+ region +'_'+ source_name + '_LC_850_phasefold'+filenameoption +'.pdf'
    plt.figure()
    plt.errorbar([x/(1/detected_freq) for x in Ydates],pfluxes, yerr=noises, fmt='o',color='g',ecolor='black',capsize=errsize,markersize=symsize)
    #plt.errorbar(Ydates[-2]/(1/detected_freq),pfluxes[-2], yerr=noises[-2], fmt='o',color='blue',ecolor='black',capsize=3)
    
    plt.errorbar([(x-Obsphases*(1/detected_freq))/(1/detected_freq) for x in Ydates],pfluxes, yerr=noises,
     fmt='o',color='g',ecolor='black',capsize=errsize,markersize=symsize)
     
    for eachObsphase in range(1,Obsphases):
        plt.errorbar([(x+eachObsphase*(1/detected_freq))/(1/detected_freq) for x in Ydates],pfluxes, yerr=noises,
         fmt='o',color='g',ecolor='black',capsize=errsize,markersize=symsize)
        plt.errorbar([(x-eachObsphase*(1/detected_freq))/(1/detected_freq) for x in Ydates],pfluxes, yerr=noises,
         fmt='o',color='g',ecolor='black',capsize=errsize,markersize=symsize)
    if highlight:
        plt.errorbar((Ydates[-1]-Obsphases*(1/detected_freq))/(1/detected_freq),pfluxes[-1], yerr=noises[-1],
        fmt='o',color='blue',ecolor='black',capsize=3)
        for eachObsphase in range(1,Obsphases):
            plt.errorbar((Ydates[-1]+eachObsphase*(1/detected_freq))/(1/detected_freq),pfluxes[-1], yerr=noises[-1],
             fmt='o',color='blue',ecolor='black',capsize=3)
            plt.errorbar((Ydates[-1]-eachObsphase*(1/detected_freq))/(1/detected_freq),pfluxes[-1], yerr=noises[-1],
             fmt='o',color='blue',ecolor='black',capsize=3)

    plt.xlim([-0.2, 1.2])
    plt.plot([x/(1/detected_freq) for x in grid], fitted_sine, color= 'Darkturquoise')
    plt.xlabel('Phase')
    plt.ylabel('Peak Flux [Jy/beam]')
    plt.title('Phase Diagram of 850 μm ('+ source_name + ' in '+ region +')', color='black')
    plt.annotate("~" + str(int(np.round(detected_period)))+ " days",xy=[0.15,0.77], xycoords="figure fraction") 
    if y_reverse:
        plt.gca().invert_yaxis()
    plt.savefig(resultpwd)
    plt.close()
    return(resultpwd)

def TCPeriodogram(source_dict, region):
    import astropy.io.ascii as astread
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.stats import LombScargle as LS
    from datetime import date
    
    i = 0
    xx = 0
    xxx = len(snames)
    det_lev = 0.0054*2
    out = open(datadir+"/Results/result.txt","w")
    print(len(snames))
    print(len(snames),file = out)
    for k in snames[0:-1]:
        pfluxes = dic[k]['peakfluxes']
        if i == 0:
            strdates = dic[k]['dates_reg']
            #print(strdates)
            ddates = []
            for i in range(len(strdates)):
                ddates.append(date(int(strdates[i][0:4]), int(strdates[i][4:6]),int(strdates[i][6:8]))
                              -date(int(strdates[i-1][0:4]), int(strdates[i-1][4:6]),int(strdates[i-1][6:8])))
            ddates[0] = (date(int(strdates[0][0:4]), int(strdates[0][4:6]),int(strdates[0][6:8]))
                             -date(int(strdates[0][0:4]), int(strdates[0][4:6]),int(strdates[0][6:8])))
            rdates = []
            for j in range(len(ddates)):
                rdates.append(ddates[j].days)
        
            for k0 in range(1,len(rdates)):
                rdates[k0] = (float(rdates[k0])+float(rdates[k0-1]))
                
            for l in range(len(rdates)):
                rdates[l] = rdates[l]/365.24
            i = i + 1
            Ydates = rdates
            #print(Ydates)
        s2n_check = 1
        for iii in range(len(pfluxes)):
            if (pfluxes[iii] < noises[iii]*3):
                s2n_check = 0
        if not s2n_check:
            xxx = xxx -1
            continue
        
        ls_pf = LS(Ydates,  pfluxes, noises)
        freq, power = ls_pf.autopower()
        fap = ls_pf.false_alarm_probability(power)
        dic[k].update({'powers':[],'frequencies':[],'faps':[]})
        dic[k]['powers'].append([a for a in power])
        dic[k]['frequencies'].append([a for a in freq])
        dic[k]['faps'].append([a for a in fap])
        if any(fap[np.where(freq < 5)] < det_lev):
            print(dic[k]['class'], dic[k]["distance"])
            print(dic[k]['class'], dic[k]["distance"], file = out)
            print(k)
            print(k, file = out)
            #print(freq[np.where(fap[np.where(freq < 5)] < det_lev)])
            #print(freq[np.where(fap[np.where(freq < 5)] < det_lev)], file = out)
            '''
            detected_freq = freq[np.where(power == power.max())]
            def func(x, a, b, c):
                return a*np.sin(b+2*np.pi/detected_freq*x)+c
            grid = np.linspace(0, 4, 1000)
            print(freq[np.where(power == power.max())])
            
            fitfunc = lambda p, x: p[0]*np.cos(2*np.sin(2*np.pi))       
            fit_result =opt.curve_fit(func, pfluxes, Ydates,np.array([0.1, 0, 1]),noises)
            xi_amp = fit_result[0][0]
            xi_xd = fit_result[0][1]
            xi_yd = fit_result[0][2]
            fitted_sine = xi_amp*np.sin(xi_xd+2*np.pi/detected_freq*grid)+xi_yd
            '''   
            cut_power = power[np.where(freq < 5)]
            detected_freq = freq[np.where(cut_power == cut_power.max())]
            detected_period = int(365.24/freq[np.where(power == power.max())])
            grid = np.linspace(-1, 4, 1000)
            fitted_sine = ls_pf.model(grid,float(detected_freq))
            plt.figure()
            plt.errorbar(Ydates,pfluxes,yerr=noises,fmt='o',color='g',ecolor='black',capsize=3)
            plt.xlim([-0.2,3.5])
            plt.plot(grid, fitted_sine, color = 'khaki')
            plt.xlabel('Relative Dates [Year]')
            plt.ylabel('Peak Flux [Jy/beam]')
            plt.title('Light Curve of ' + wl + 'μm ('+ k + ' in '+ region[regionind] +')')
            #if not(detected_freq == freq.min()):
            plt.annotate("~" + str(detected_period).format("4f")+ " days",xy=[0.15,0.77], xycoords="figure fraction") 
            plt.savefig(datadir+'/Results/'+ region[regionind] +'_'+ k + '_LC_' + wl + '.pdf')
            plt.close()
            plt.figure()
            p_thr = ls_pf.false_alarm_level(det_lev)
            p_tough_thr = ls_pf.false_alarm_level(0.001)
            plt.plot(freq, power,color='g')
            plt.plot([-1,6],[p_thr, p_thr], color='black')
            plt.plot([-1,6],[p_tough_thr, p_tough_thr], color='gray')
            #if not(detected_freq == freq.min()):
            plt.text(freq[np.where(power == power.max())], power.max()+0.05, "~" +str(detected_period) + ' days')
            plt.xlabel("Cycles per year")
            plt.ylabel("Power")
            plt.title("Power Spectrum of " + wl + 'μm (' + k +' in '+ region[regionind] +')')
            plt.xlim([0,5])
            plt.ylim([0,1])
            plt.savefig(datadir+'/Results/'+ region[regionind] +'_'+ k + '_PS_' + wl + '.pdf')
            plt.close()
            xx = xx +1
    print(xxx)
    print(xxx, file = out)
    print(xx)
    print(xx, file = out)
    out.close()
    np.save(datadir+'/'+region[regionind]+'_peak_track_withnoise_withPS.pkl.npy',dic)s


def TCCalNoiseadd(noises, pfluxes):
    import numpy as np
    
    caladdnoise = np.sum([noises,[0.02*x for x in pfluxes]],axis = 0)
    return(caladdnoise)
    
    
def TCShuffling(pfluxes, Ydates, noises, N_iter=1000):
    import numpy as np
    from astropy.stats import LombScargle as LS
    import random

    n = len(Ydates)
    #print(pfluxes)
    result = []
    Power_maximums = []
    Frequency_maximums = []
    a = tuple(pfluxes)
    print("\n")
    for i in range(0,N_iter):
        random.shuffle(pfluxes)
        ls_pf = LS(pfluxes, Ydates, noises)
        Freq, Power = ls_pf.autopower()
        Power_maximums.append(Power.max())
        Frequency_maximums.append(Freq[np.where(Power == Power.max())])
        #print(pfluxes)
        result.append([])
        for j in range(0, n):
            result[i].append(pfluxes[j])
        #print(result)
        if i%1000 == 0:
            print(i)
    return(Power_maximums, Frequency_maximums)
    

def PlanckF(x,T):
    import numpy as np
    h =6.626e-34    #[m^2 kg s-1]
    c =2.998e+8      #[m s-1] 
    k =1.381e-23    #[m^2 kg s-2 K]   
    
    return(((2*h*(c**2))/(x**5))/(np.exp((h*c)/(x*k*T))-1))

def SEDflux(T):
    import numpy as np
    import scipy.integrate as integrate
    
    bands = {'U':[320e-9,400e-9],'B':[400e-9,500e-9],'V':[500e-9,700e-9],'R':[550e-9,800e-9],'I':[700e-9,900e-9]}
    filters = list(bands.keys())
    #x1 =10^(-10)
    #x2 =10^(-3.1)  ; Integration near to zero to near to infinite
    results_bands = []
    for band in filters:
        x1 = bands[band][0]
        x2 = bands[band][1]
        flux = integrate.quad(lambda x: func(x,T), x1, x2)[0]
        results_bands.append(flux)
        print('Intensity in ' + band + ' band is '+ str(flux))
    results=np.array(results_bands)
    maximum = np.where(results == results.max())

    print('maximum among the UBVRI filters is '+filters[maximum[0][0]])

    x1 =4*10^(-7)
    x2 =7*10^(-7)
    result_quad = integrate.quad(lambda x: func(x, T), x1, x2)
    print('scipy.integrate.quad Result = '+str(result_quad[0]))
    
    
    
    
