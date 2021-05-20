def JCMTTransient(region,mod_date,wl='850'):
    import os
    import numpy as np
    import TCTriggerFunctions as TF
    
    import importlib
    importlib.reload(TF)


    
    #region = ["IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", "SERPS"]
    #FieldName is one of: "IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", or "SERPS"

    datadir = "/mnt/Secdrive/TSmain/"+region
    if wl == '450':
        datadir = "/mnt/Secdrive/TSmain/"+region + '/450'
    output_dir = datadir + '/Results_' + mod_date +'/'
    os.system('mkdir ' + datadir + '/' + 'Results_'+mod_date)
    dic =np.load(datadir+'/'+region +'_peak_noise_YSOcom_'+mod_date+'.pkl.npy').item()
    noises = dic['noises']
    option = 'ACFtest'
    snames = list(dic.keys())
        
    return(snames[:-1], dic)
    
class JCMTTScall:

    def __init__(self, rgs, index, ID=0):
        import Table_Reader as TR
        import numpy as np
        input_dir = "/media/lyhee/Secdrive/TSmain/TSmain/"
        self.data = TR.ReadTable(input_dir+'Data/'+rgs+'_source_data.dat')
        self.JDs = self.data.JDs
        self.pfluxes = self.data.pfluxes[index]
        self.sd_fids = self.data.sd_fids[index]
        if ID:
            self.ID = ID
        if not ID:
            self.ID = self.data.IDs[index]
        noises = np.full(len(self.pfluxes),self.sd_fids)
        self.mags, self.mnoises = Mag2flux("",self.pfluxes,noises,unit="Jy",reverse=1)
        
    def matchWISE(self):
        self.Windex=0
        
        if self.ID == 'JCMTPP_J034357.0+320305':
            self.Windex=6095
        if self.ID == 'JCMTPP_J032910.4+311331':
            self.Windex=5953
        if self.ID == 'JCMTPP_J032903.8+311449':
            self.Windex=5933
        if self.ID == 'JCMTPP_J032911.1+311828':
            self.Windex=5958
        #if self.ID == 'JCMTPP_J054607.2-001332':
        #    self.Windex=3161
        if self.ID == 'JCMTPP_J054603.6-001447':
            self.Windex=3159
            
        if self.ID == 'JCMTPP_J054613.2-000602':
            self.Windex=3180
            
        #if self.ID == 'JCMTPP_J053522.4-050111':
        #    self.Windex=2437
        if self.ID == 'JCMTPP_J053529.8-045944':
            self.Windex=2458
        
        if self.ID == 'JCMTPP_J182949.8+011520':
            self.Windex=6375
        if self.ID == 'JCMTPP_J182948.2+011644':
            self.Windex=6369
        
        if self.ID == 'JCMTPP_J182952.0+011550':
            self.Windex=6384
            
        return(self.Windex)
        

def JCMTTRansientSource(k,region,wl,mod_date,timeunit="day",test=0):
    import importlib
    import os
    import numpy as np
    import TCTriggerFunctions as TF
    import Table_Reader as TR
    importlib.reload(TR)
    importlib.reload(TF)
    
    #region = ["IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", "SERPS"]
    #FieldName is one of: "IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", or "SERPS"
    
    datadir = "/mnt/Secdrive/TSmain/"+region
    
    if test:
        data = TR.ReadTable(datadir+'/Data/'+region+'_source_info_'+mod_date[-4:]+'.dat')

        if wl == '450':
            data = TR.ReadTable(datadir+'/Data/450'+region+'_source_info_'+mod_date[-4:]+'.dat')

        k = TF.TCCelebrities(k)
        if k[:4] != 'JCMT':
            k = TF.TCCelebrities(k)
        i = data.FindIndex(ID=k)
        
        return(data.JDs, data.pfluxes[i], data.sd_map)

    
    output_dir = datadir + '/Results_' + mod_date +'/'
    os.system('mkdir ' + datadir + '/' + 'Results_'+mod_date)
    if wl == '850':
        dic =np.load(datadir+'/'+region+'_peak_noise_YSOcom_'+mod_date+'.pkl.npy').item()
    if wl == '450':
        dic = np.load(datadir+'/450/'+region+'_peak_noise_YSOcom_'+mod_date+'.pkl.npy').item()
    noises = dic['noises']
    option = 'ACFtest'
    snames = list(dic.keys())
    
    k = TF.TCCelebrities(k)
    if k[:4] != 'JCMT':
        k = TF.TCCelebrities(k)
    pfluxes = dic[k]['peakfluxes']
    noises = dic['noises']
    if wl == '850':
        noises = np.sum([noises, [0.02*x for x in pfluxes]], axis=0)
    if wl == '450':
        noises = np.sum([noises, [0.05*x for x in pfluxes]], axis=0)
    JDs = dic[k]['dates']

    dates = [(x - JDs[0]) for x in JDs]
    if timeunit == "year":
        dates = [(x - JDs[0])/365.24 for x in JDs]
    Start_date = dic[k]['dates_reg'][0]
    Last_date = dic[k]['dates_reg'][-1]
    
    
        
    return(np.asarray(JDs), np.asarray(pfluxes), np.asarray(noises))
        
def EC53IR(wl,imeunit="day",addLiverpool=0,timeunit="day"):
    import os
    import numpy as np
     #===========================IR Data Call=====================================
    #FieldName is one of: "IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", or "SERPS"
    datadir = "/mnt/Secdrive/TSmain/SERPM"
    k = 'EC53'
    #option = "_Yoo450"
    option = ""
    output_dir = datadir + '/EC53/EC53'
    input_dir = output_dir
    wlsplit = wl.split("_")
    if wlsplit[0] == "Hodapp":
        JDs,mags = np.loadtxt(input_dir +"/EC53_"+wl+".dat", unpack=True, dtype='str')
        mags = [float(e) for e in mags]
        JDs = [float(e) for e in JDs]
        noises = [0.05 for e in JDs]
    if wlsplit[0] == "UKIRT":
        JDs,mags,noises = np.loadtxt(input_dir +"/EC53_"+wl+".dat", unpack=True, dtype='str')
        #[JDs,regdates,mags,noises,temp] = np.loadtxt(input_dir +"/EC53_"+wl+".dat", unpack=True, dtype='str')
        #JDs_new, mags_new, noises_new = np.loadtxt(input_dir+"/ec53_7.mean",unpack=True,dtype='str')
        #regdates = [e[1:] for e in regdates]
        mags = [float(e) for e in mags]
        JDs = [float(e) for e in JDs]
        if wlsplit[0] == 'J':
            JDs_begin = JDs[0]  #2456942.789833
        JDs_begin = 2456942.789833
        noises = [float(e) for e in noises]

        #Ydates = SAn.TCConvertDates(regdates)
        Ydates = [(x - JDs_begin) for x in JDs]

    if wlsplit[0] == "Liverpool":
        JDs,mags,noises,NPTS = np.loadtxt(input_dir +"/EC53_Liverpool_H.dat", unpack=True, dtype='float')
        JDs = [float(e)+2450000.5 for e in JDs]
        mags = [float(e)+22.93 for e in mags]
        noises = [float(e) for e in noises]
        
    if wlsplit[0] == "IRIS":
        JDs,mags,noises = np.loadtxt(input_dir +"/EC53_"+wl+".dat", unpack=True, dtype='str')
        mags = [float(e) for e in mags]
        JDs = [float(e) for e in JDs]
        noises = [float(e) for e in noises]
        
    if wlsplit[0] == "Combined":
        JDs,mags,noises = np.loadtxt(input_dir +"/EC53_"+wl+".dat", unpack=True, dtype='str')
        mags = [float(e) for e in mags]
        JDs = [float(e) for e in JDs]
        noises = [float(e) for e in noises]
    
    if addLiverpool:
        if wlsplit[1] == "H":
            JDs_LP,mags_LP,noises_LP = np.loadtxt(input_dir +"/EC53_Liverpool_H.dat", unpack=True, dtype='float')
            JDs_begin = 2456942.789833
            Ydates_LP = [float(e)+2450000.5-JDs_begin for e in JDs_LP]
            JDs_LP = [float(e)+2450000.5 for e in JDs_LP]
            mags_LP = [float(e)+0.15 for e in mags_LP]
            noises_LP = [float(e) for e in noises_LP]
            for everydate in JDs_LP:
                JDs.append(everydate)
            for everymag in mags_LP:
                mags.append(everymag)
            for everynoise in noises_LP:
                noises.append(everynoise)
            mags_copy = tuple(mags)
            noises_copy = tuple(noises)
            for eachindex in range(0,len(JDs)):
                sortedindex = sorted(JDs).index(JDs[eachindex])
                mags[sortedindex] = mags_copy[eachindex]
                noises[sortedindex] = noises_copy[eachindex]
            JDs.sort()
            Ydates = [(x - JDs_begin) for x in JDs]
    if wlsplit[0] == "WISE":
        JDs,mags,noises,band,dist,ra,dec = np.loadtxt(input_dir +"/EC53_WISE.asc", unpack=True, dtype='str')
        mags = np.array([float(e) for e in mags])
        JDs = np.array([float(e) for e in JDs])
        noises = np.array([float(e) for e in noises])
        if wlsplit[2] == "all":
            JDs,mags,noises,band,dist,ra,dec = np.loadtxt(input_dir +"/EC53_WISE_all.asc", unpack=True, dtype='str')
        if wlsplit[1] == "W1":
            JDs = JDs[np.where(band == "W1")]
            mags = mags[np.where(band == "W1")]
            noises = noises[np.where(band == "W1")]
        if wlsplit[1] == "W2":
            JDs = JDs[np.where(band == "W2")]
            mags = mags[np.where(band == "W2")]
            noises = noises[np.where(band == "W2")]
        if wlsplit[1] == "W3":
            JDs = JDs[np.where(band == "W3")]
            mags = mags[np.where(band == "W3")]
            noises = noises[np.where(band == "W3")]
        if wlsplit[1] == "W4":
            JDs = JDs[np.where(band == "W4")]
            mags = mags[np.where(band == "W4")]
            noises = noises[np.where(band == "W4")]
    if wl == 'Yoo_K':
        noises = [0.05 for x in JDs]
    if timeunit=="year":
        Ydates = [x/365.24 for x in Ydates]

    return(np.array(JDs),np.array(mags),np.array(noises))

class WISEcall:
    def  __init__(self, index, wv,start=1):
        self.inputdir= "/mnt/Secdrive/TSmain/TSmain/WISE/"
        self.filename= "_cavg"
        self.wv = '_'+ wv
        if self.wv == '_W2':
            self.wv = ''
        self.index = index
        self.inputfile = self.inputdir+str(self.index)+self.filename+self.wv+'.csv'
        order, JDs, mags, noises = np.loadtxt(self.inputfile, unpack=True, dtype='str',delimiter=',')
        self.JDs = [float(x)+2400000.5 for x in JDs[1:]]
        self.mags = [float(x) for x in mags[start:]]
        self.noises = [float(x) for x in noises[start:]]
        #JDs = JDs[np.where(band == "W2")]

        #self.mags = mags[np.where(band == "W2")]
        #self.noises = noises[np.where(band == "W2")]


def Mag2flux(wl,mags,noises,unit="Jy",reverse=0):
    import numpy as np
    wlsplit = wl.split("_")

    if wlsplit[0] == "UKIRT":
        if wlsplit[1] == "J":
            Z = 1600
        if wlsplit[1] == "H":
            Z = 1020
        if wlsplit[1] == "K":
            Z = 657
    if wlsplit[0] == "Liverpool":
            Z = 1035.9
    if wlsplit[0] == "WISE":
        if wlsplit[1] == "W1":
            Z = 306.681
        if wlsplit[1] == "W2":
            Z = 170.663
        if wlsplit[1] == "W3":
            Z = 29.0448
        if wlsplit[1] == "W4":
            Z = 8.2839
        #if wlsplit[0] == "Hodapp":
        
    if reverse:
        fnoises = [x*(2.5)/y/np.log(10) for x,y in zip(noises,mags)]
        fluxes = [(-2.5)*np.log10(x/mags[0]) for x in mags]
        #if wlsplit[0] != "850":
        #    fluxes = [(-2.5)*np.log10(x/Z) for x in mags]
        return(np.array(fluxes),np.array(fnoises))
        
    fluxes = [10**(x/(-2.5))*Z for x in mags]
    fnoises = [x*(1/2.5)*np.log(10)*y for x,y in zip(fluxes,noises)]

    if unit == "mJy":
        fluxes = [x*1000 for x in fluxes]
        fnoises = [x*1000 for x in fnoises]
        
    return(np.array(fluxes),np.array(fnoises))


def LCexport(k,regionind,dates,fluxes,noises,mod_date, wl="850"):
    region = ["IC348", "NGC1333", "NGC2024", "NGC2068", "OMC23", "OPHCORE", "SERPM", "SERPS"]
    datadir = "/mnt/Secdrive/TSmain/"+region[regionind]
    outputdir = datadir+"/Results_" + mod_date + "/" + k +"_850_"+mod_date+".txt"
    out = open(outputdir,"w")
    print("#   " + k +" by JCMT Transient Survey" ,file = out)
    print("#   " + "850 micron (SCUBA2)",file = out)
    print("#   " + "Updated on " + mod_date, file =out)
    print("   MJD    Flux[Jy/beam]    Noise[Jy/beam]",file=out)
    for i in range(len(dates)):
        print(str(dates[i]) +" "+ str(fluxes[i]) +" "+ str(noises[i]),file=out)
    out.close()
    
    return(outputdir)
        
def JD2date(mjd, fmt='mjd'):
    import julian
#    dates = []
    dt = julian.from_jd(mjd, fmt='mjd')
    if fmt == 'jd':
        dt = julian.from_jd(mjd, fmt='jd')
        
    date = str('{:04d}').format(dt.year) + str('{:02d}').format(dt.month) +str('{:02d}').format(dt.day)
    return(date)

def Combinedata(JD1, data1, noise1, JD2, data2, noise2):
    import numpy as np
    JD1, data1, noise1 = list(JD1), list(data1), list(noise1)
    JD2, data2, noise2 = list(JD2), list(data2), list(noise2)
    for everydate in JD2:
        JD1.append(everydate)
    for everydata in data2:
        data1.append(everydata)
    for everynoise in noise2:
        noise1.append(everynoise)
    data1_copy = tuple(data1)
    noise1_copy = tuple(noise1)
    for eachindex in range(0,len(JD1)):
        sortedindex = sorted(JD1).index(JD1[eachindex])
        data1[sortedindex] = data1_copy[eachindex]
        noise1[sortedindex] = noise1_copy[eachindex]
    JD1.sort()
    return(np.array(JD1), np.array(data1), np.array(noise1))
        
#def LCHyperthetic(regionind, 

import numpy as np
class EC53call:
    import LCCall as LCC
    import importlib
    importlib.reload(LCC)
    k = "EC53"
    mod_date='191021'
    dates_K, mags_K, noises_K = LCC.EC53IR("UKIRT_K",imeunit="day",addLiverpool=0,timeunit="day")
    dates_IK, mags_IK, noises_IK = LCC.EC53IR("IRIS_Ks",imeunit="day",addLiverpool=0,timeunit="day")
    #dates_CK, mags_CK, noises_CK = LCC.EC53IR("Combined_K",imeunit="day",addLiverpool=0,timeunit="day")
    dates_CK, mags_CK, noises_CK = LCC.Combinedata(dates_K,mags_K,noises_K,dates_IK,mags_IK,noises_IK)
    
    dates_H, mags_H, noises_H = LCC.EC53IR("UKIRT_H",imeunit="day",addLiverpool=0,timeunit="day")
    dates_IH, mags_IH, noises_IH = LCC.EC53IR("IRIS_Hn",imeunit="day",addLiverpool=0,timeunit="day")
    dates_LH, mags_LH, noises_LH = LCC.EC53IR("Liverpool_H",imeunit="day",timeunit="day")
    #dates_CH, mags_CH, noises_CH = LCC.EC53IR("Combined_H",imeunit="day",addLiverpool=0,timeunit="day")
    dates_IUH, mags_IUH, noises_IUH = LCC.Combinedata(dates_H,mags_H,noises_H,dates_IH,mags_IH,noises_IH)
    dates_CH, mags_CH, noises_CH = LCC.Combinedata(dates_IUH,mags_IUH,noises_IUH,dates_LH,mags_LH,noises_LH)
    
    dates_J, mags_J, noises_J = LCC.EC53IR("UKIRT_J",imeunit="day",addLiverpool=0,timeunit="day")
    dates_IJ, mags_IJ, noises_IJ = LCC.EC53IR("IRIS_Jn",imeunit="day",addLiverpool=0,timeunit="day")
    dates_CJ, mags_CJ, noises_CJ = LCC.EC53IR("Combined_J",imeunit="day",addLiverpool=0,timeunit="day")
    dates_CJ, mags_CJ, noises_CJ = LCC.Combinedata(dates_J,mags_J,noises_J,dates_IJ,mags_IJ,noises_IJ)

    dates_850, fluxes_850, noises_850 = LCC.JCMTTRansientSource(k,"SERPM","850",'0410',test=1)
    dates_450, fluxes_450, noises_450 = LCC.JCMTTRansientSource('JCMTPP_J182951.3+011638',"SERPM","450","0410",test=1)
    JDs_W1, mags_W1, noises_W1 = LCC.EC53IR("WISE_W1_avg")
    JDs_W2, mags_W2, noises_W2 = LCC.EC53IR("WISE_W2_avg")
    JDs_W3, mags_W3, noises_W3 = LCC.EC53IR("WISE_W3_avg")
    JDs_W4, mags_W4, noises_W4 = LCC.EC53IR("WISE_W4_avg")
    dates_HK, mags_HK, noises_HK = LCC.EC53IR("Hodapp_K")
    mags_850,mnoises_850 = LCC.Mag2flux("850",fluxes_850,noises_850, reverse=1)
    mags_450,mnoises_450 = LCC.Mag2flux("450",fluxes_450,noises_450, reverse=1)
    input_dir = '/mnt/Secdrive/TSmain/SERPM/EC53/EC53/'
    JDs,Js,Std_Js,Hs,Std_Hs,Ks,Std_Ks,J_Hs,J_H_errs,H_Ks,H_K_errs, = np.loadtxt(input_dir +"/EC53_UKIRT_JHK_combined_fin.dat", unpack=True, dtype='float')
    
