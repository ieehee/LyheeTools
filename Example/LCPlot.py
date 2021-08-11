
def LCplot(Ydates, pfluxes, noises,mod_date,region,k,wl='850',mag=0,option='',fitting=0,string=0,save=0,newplot=0
          , color='green',ecolor = 'black',stringcolor='green',O_F=0,figureindex=1,markersize=5,capsize=3,alpha=1,label="",
          legend=0,title=0,xlabel='MJD'):
    import matplotlib.pyplot as plt
    import numpy as np
    
    if newplot:
        plt.close()
        plt.figure(figureindex)
    #pfluxes = [x*1000 for x in pfluxes]
    plt.errorbar(Ydates,pfluxes,yerr=noises,fmt='o',markersize=markersize,color=color,ecolor=ecolor,
                label = label,capsize=capsize,alpha=alpha)
    if string:
        plt.plot(Ydates,pfluxes,color=stringcolor)
    #plt.xlabel('Relateive Dates [days]')
    plt.xlabel(xlabel,fontsize=13)
    plt.ylabel('Flux [mJy]',fontsize=13)
    if wl == "850":
        plt.ylabel('Flux [Jy/beam]')
#    plt.ylabel("Flux $\cross$  [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]")
        
    if mag:
        plt.ylabel('Magnitude [mag]',fontsize=15)
        plt.gca().invert_yaxis()
    if title:
        plt.title('Light Curve of ' + k + ' in ' + wl)
    datadir = "/mnt/Secdrive/TSmain/"+region
    output_dir = datadir+"/Results_" + mod_date +"/"
    if k == "EC53":
        output_dir = output_dir + "EC53/"
    output = output_dir+ region +'_'+ k + '_LC_' + wl +option+ '.pdf'
    if O_F:
        output = output_dir+ region +'_'+ k + '_OF_' + wl +option+ '.pdf'
        plt.plot([-100,np.array(Ydates).max()+100],[0,0],color=stringcolor)
    if legend ==1:
        plt.legend(loc = "upper left")
    if legend ==2:
        plt.legend(loc = "upper right")
    if save:
        plt.savefig(output)

    return(output)
    

def JDUTPlot(ax, MJD=1, fs=13,nolabel=0,spacing="      ",noedge=[0,0]):
    # Overplot the UT notation over the JD lightcurve plot
    import matplotlib.pyplot as plt
    import numpy as np
    
    MJDcoverage = np.array([x[0] for x in ax.dataLim.get_points()])
    #MJDrange = MJDcoverage[1]-MJDcoverage[0]
    
    UTref=1995
    if MJD:
        MJDref = 49718 # 1995
    else:
        MJDref = 2449718.5
        
    MJD_startline = MJDcoverage[0]-MJDref # JD_begin - JD(1995 yr)
    MJD_endline = MJDcoverage[1]-MJDref # JD_end - JD(1995 yr )
    Yr_fromref = int(MJD_startline/365)
    
    UT_startline = Yr_fromref +UTref+1 # JD_begin[yr] + UTref(1995) +1
    UT_endline = int(MJD_endline/365)+UTref+1
    
    MJD_starter = int(MJD_startline/365+1)*365+MJDref
    
    #MJD_starter = (UT_startline-int(UT_startline))*365
    
    #=UT_startline - UT_endline
    #print("MJD_startline: {}".format(MJD_startline))
    #print("MJD_starter: {}".format(MJD_starter))
    UT_gridpos = np.arange(MJD_starter, MJD_endline+MJDref,365)
    #print(UT_gridpos)
    UT_grid = np.arange(UT_startline, UT_endline,1)
    #print(UT_grid)
    UT_labels = [spacing+str(x) for x in UT_grid]#+[' ']
    if noedge[0]:
        UT_labels[0]=' '
    if noedge[1]:
        UT_labels[-1]=' '

    ax2 = ax.twiny()
    #setxticks = ax.get_xticks()
    ax2.set_xlim(MJDcoverage[0], MJDcoverage[1])
    ax2.set_xticks(UT_gridpos)
    ax2.grid(True,color='grey',alpha=0.3)
    if nolabel:
        ax2.set_xticklabels([' ' for x in UT_labels])
    else:
        ax2.set_xlabel('UT',fontsize=fs+1)
        ax2.set_xticklabels(UT_labels,fontsize=fs,horizontalalignment='left')
    
    return(ax2)
    

def ACFplot(r_dates, r_fluxes, r_noises,ACF, mod_date, wl,k,region, interval,revolving,save=0,option='',title=0,mag=0):
    import matplotlib.pyplot as plt
    import numpy as np
    
    r_dates_plot = r_dates[0:int(len(r_dates)/2)]
    fig, ax = plt.subplots(2,1,figsize=(5,8))
    start_ind = int(300/interval)
    j = start_ind + np.where(ACF[start_ind:int(len(r_dates)/2)] == ACF[start_ind:int(len(r_dates)/2)].max())[0]
    print(j)
    j = j[0]
    #j = 16
    i = len(r_dates)-j
    if revolving:
        ax[0].errorbar(r_dates,list(r_fluxes)[i:]+list(r_fluxes)[:i],yerr=list(r_noises)[i:]+list(r_noises)[:i], fmt='o',     color='orange')
    if not revolving:
        ax[0].errorbar(r_dates[j:],list(r_fluxes)[:i],yerr=list(r_noises)[:i], fmt='o', color='orange')
    ax[0].errorbar(r_dates, r_fluxes, yerr=r_noises, fmt='o', color='firebrick')
    if mag:
        ax[0].invert_yaxis()
    ax[0].set(xlabel = 'Relative Dates [days]',ylabel = 'Flux [mJy]')
    if wl == "850":
        ax[0].set(xlabel = 'Relative Dates [days]',ylabel = 'Peak Flux [Jy/beam]')
    
    if mag:
        ax[0].set(xlabel = 'Relative Dates [days]',ylabel = 'Magnitudes [mag]')
    if title:
        #ax[0].set_title('Light Curve of ' + wl + ' ('+ k + ' in '+ region +')')#plt.close()
        ax[0].annotate("Resampled \n"+title,xy=[0.2,0.83], xycoords="figure fraction")
    ax[1].scatter(r_dates_plot,ACF)
    ax[1].plot([j*interval,j*interval],[ACF.min(),ACF.max()],color="orange")
    ax[1].annotate("Peak at {0:3.0f}$\pm${1:2.1f} days".format(j*interval,interval/2),xy=[0.40,0.35], xycoords="figure fraction")
    ax[1].set(xlabel = 'Shifted Days [days]', ylabel = 'ACF')
    datadir = "/mnt/Secdrive/TSmain/"+region
    output_dir = datadir+"/Results_" + mod_date + "/Period/"
    output = output_dir+ region +'_'+ k + '_ACF_' + wl +option+"_"+str(interval)+ 'int_check.pdf'
    
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    
    if save:
        plt.savefig(output)
        
    return(output)
    
def LSPplot(x, power,xunit="period"):
    import matplotlib.pyplot as plt
    import numpy as np
    
    if xunit=="period":
    
        plt.xlabel("Period [day]")
    plt.ylabel
    return(output)
    
def LCPhaseFold(source_name,region,pfluxes,noises,grid,fitted_sine,Ydates,detected_freq,datadir,filenameoption,highlight=0,y_reverse=0,errsize=3,symsize=5):
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
    plt.annotate("~" + str(int(detected_period))+ " days",xy=[0.15,0.77], xycoords="figure fraction") 
    if y_reverse:
        plt.gca().invert_yaxis()
    plt.savefig(resultpwd)
    plt.close()
    return(resultpwd)

