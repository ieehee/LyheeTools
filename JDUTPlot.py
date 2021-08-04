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
    
    UT_gridpos = np.arange(MJD_starter, MJD_endline+MJDref,365)
    #print(UT_gridpos)
    UT_grid = np.arange(UT_startline, UT_endline,1)
    #print(UT_grid)
    UT_labels = [spacing+str(x) for x in UT_grid]
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
