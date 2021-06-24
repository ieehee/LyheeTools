import numpy as np
import LCAnalyses as LCA

class ReadTable:
    def __init__(self,inputfile,raw=0,datecut=0,recal=0):
        # Call the lightcurve of the sources in a region
        # input: path to the source_info table.
        # Requires the metadata file matches to the source_info table, in the same directory
        self.inputfile = inputfile
        self.a = np.loadtxt(self.inputfile,dtype='str')
        self.metafile = inputfile.split('_')[0] + '_meta.dat'
        self.b = np.loadtxt(self.metafile, dtype='str')
        self.Indices = [int(x[0]) for x in self.a]
        self.UTs = [x[2] for x in self.b]
        self.JDs = [float(x[4]) for x in self.b]
        if recal:
            self.JDs = [float(x[3]) for x in self.b]
        self.IDs = [x[1] for x in self.a]
        self.RAs = [np.float64(x[2]) for x in self.a]
        self.DECs = [np.float64(x[3]) for x in self.a]
        self.proto_dists = [np.float16(x[4]) for x in self.a]
        self.disk_dists = [np.float16(x[5]) for x in self.a]
        self.mean_pfluxes= [np.float64(x[6]) for x in self.a]
        self.sd_map = [np.float64(x[7]) for x in self.b]
        #self.sd_map = [x[8] for x in self.b]
        self.sd_pfluxes = [np.float64(x[7]) for x in self.a]
        self.sd_fids = [np.float64(x[8]) for x in self.a]
        #self.slopes = [np.float64(x[10]) for x in self.a]
        #self.dslopes = [np.float64(x[11]) for x in self.a]
        self.pfluxes = [np.float64(x[15:15+len(self.b)]) for x in self.a]
        self.deviations = [np.float64(x[14+len(self.b):]) for x in self.a]
        #print(len(self.deviations[0]))
        
        # Below _if_ blocks are to exclude the problematic epochs in NGC1333 and OMC23
        if inputfile.split("_")[0].split("/")[-1] == "NGC1333" and not raw:
            #print('Occured!')
            for index in [2,15,18]: # 3, 17, 21th
                #print(self.UTs[index])
                self.JDs = np.delete(self.JDs, [index])
                self.UTs = np.delete(self.UTs, [index])
                for i in range(0,len(self.pfluxes)):
                    self.pfluxes[i] = np.delete(self.pfluxes[i], [index])
                    self.deviations[i] = np.delete(self.deviations[i], [index])
                    
                #print(len(self.deviations[0]))
            self.mean_pfluxes = np.array([np.mean(x) for x in self.pfluxes])
            self.sd_pfluxes = [np.std(x) for x in self.pfluxes]
            self.sd_fids = np.sqrt(0.014**2 + (0.02*self.mean_pfluxes)**2)
            
        elif inputfile.split("_")[0].split("/")[-1] == "OMC23" and not raw:
            #print("Occured!")
            index = 18
            self.JDs = np.delete(self.JDs, [index])
            for i in range(0,len(self.pfluxes)):
                self.pfluxes[i] = np.delete(self.pfluxes[i], [index])
                self.deviations[i] = np.delete(self.deviations[i], [index])
            #print(len(self.deviations[0]))
            self.mean_pfluxes = np.array([np.mean(x) for x in self.pfluxes])
            self.sd_pfluxes = [np.std(x) for x in self.pfluxes]
            self.sd_fids = np.sqrt(0.014**2 + (0.02*self.mean_pfluxes)**2)
            
        # For testing the recalibrated (Steve & Colton et al. in prep) data.
        if recal:
            self.recal_factor = [np.float64(x[-1]) for x in self.b]
            for i in range(0,len(self.pfluxes)):
                for j in range(0,len(self.pfluxes[0])):
                    self.pfluxes[i][j] = self.pfluxes[i][j]*self.recal_factor[j]
        
        # Limit the date for the data
        if datecut and datecut < np.max(self.JDs):
            print(datecut)
            if datecut < np.min(self.JDs):
                print("Recheck the datecut")
                #aaa=aaa
            dates = np.array(self.JDs)
            cutter = np.where(datecut > dates)
            cutter = cutter[-1]
            #self.Indices = self.Indices[:cutter[-1]+1]
            self.JDs = self.JDs[:cutter[-1]+1]
            print(self.JDs[-1])
            #self.IDs = self.IDs[:cutter[-1]+1]
            #self.RAs = self.RAs[:cutter[-1]+1]
            #self.DECs = self.DECs[:cutter[-1]+1]
            #self.proto_dists = self.proto_dists[:cutter[-1]+1]
            #self.disk_dists = self.disk_dists[:cutter[-1]+1]

            #self.slopes = [np.float64(x[10]) for x in self.a]
            #self.dslopes = [np.float64(x[11]) for x in self.a]
            self.pfluxes = [x[:cutter[-1]+1] for x in self.pfluxes]
            self.deviations = [x[:cutter[-1]+1] for x in self.deviations]
            self.mean_pfluxes = np.array([np.mean(x) for x in self.pfluxes])
            self.sd_pfluxes = [np.std(x) for x in self.pfluxes]
            self.sd_fids = np.sqrt(0.014**2 + (0.02*self.mean_pfluxes)**2)
            
    
    def FindIndex(self, ID=0, Index =0):
        # Input ID -> return index
        # Input index -> return ID
        if Index:
            return(self.Indices.index(Index))
        if ID:
            return(self.IDs.index(ID))
    
    def Identifier(self, ID=0, Index =0):
        # Check whether a source has other known name
        if not ID:
            for i in range(0, len(self.pfluxes)):
                if self.IDs[i]=='JCMTPP_J034356.5+320050':
                    self.IDs[i] = "MMS 1"
                if self.IDs[i] =='JCMTPP_J034357.0+320305':
                    self.IDs[i] = "HH 211"
                if self.IDs[i] =='JCMTPP_J032910.4+311331':
                    self.IDs[i] = "IRAS 4A"
                if self.IDs[i] =='JCMTPP_J032903.4+311558':
                    self.IDs[i] = "VLA 3"
                if self.IDs[i] =='JCMTPP_J032903.8+311449':
                    self.IDs[i] = "West 40"
                if self.IDs[i] =='JCMTPP_J032911.1+311828':
                    self.IDs[i] = "HH 6"
                if self.IDs[i] =='JCMTPP_J054607.2-001332':
                    self.IDs[i] = "HOPS 358"
                if self.IDs[i] =='JCMTPP_J054603.6-001447':
                    self.IDs[i] = "HOPS 315"
                if self.IDs[i] =='JCMTPP_J054631.0-000232':
                    self.IDs[i] = "HOPS 373"
                if self.IDs[i] =='JCMTPP_J054647.4+000028':
                    self.IDs[i] = "HOPS 389"
                if self.IDs[i] =='JCMTPP_J054613.2-000602':
                    self.IDs[i] = "V1647 Ori"
                if self.IDs[i] =='JCMTPP_J053529.8-045944':
                    self.IDs[i] = "HOPS 383"
                if self.IDs[i] =='JCMTPP_J162626.8-242431':
                    self.IDs[i] = "VLA 1623-243"
                if self.IDs[i] =='JCMTPP_J182949.8+011520':
                    self.IDs[i] = "SMM 1"
                if self.IDs[i] =='JCMTPP_J182948.2+011644':
                    self.IDs[i] = "SH 2-68 N"
                if self.IDs[i] =='JCMTPP_J182951.2+011638':
                    self.IDs[i] = "EC 53"
                if self.IDs[i] =='JCMTPP_J182952.0+011550':
                    self.IDs[i] = "SMM 10"
                if self.IDs[i] =='JCMTPP_J183004.0-020306':
                    self.IDs[i] = "CARMA 7"
                if self.IDs[i] =='JCMTPP_J182937.8-015103':
                    self.IDs[i] = "IRAS 18270-0153"
                if self.IDs[i] =='JCMTPP_J053527.4-050929':
                    self.IDs[i] = "HOPS 370"
                if self.IDs[i] =='JCMTPP_J053522.4-050111':
                    self.IDs[i] = "HOPS 88"
                if self.IDs[i] =='JCMTPP_J183002.6-020248':
                    self.IDs[i] = "CARMA 3"

def Identifier(ID):
    if ID=='JCMTPP_J034356.5+320050':
        ID = "MMS 1"
    if ID =='JCMTPP_J034357.0+320305':
        ID = "HH 211"
    if ID =='JCMTPP_J032910.4+311331':
        ID = "IRAS 4A"
    if ID =='JCMTPP_J032903.4+311558':
        ID = "VLA 3"
    if ID =='JCMTPP_J032903.8+311449':
        ID = "West 40"
    if ID =='JCMTPP_J032911.1+311828':
        ID = "HH 6"
    if ID =='JCMTPP_J054607.2-001332':
        ID = "HOPS 358"
    if ID =='JCMTPP_J054603.6-001447':
        ID = "HOPS 315"
    if ID =='JCMTPP_J054631.0-000232':
        ID = "HOPS 373"
    if ID =='JCMTPP_J054647.4+000028':
        ID = "HOPS 389"
    if ID =='JCMTPP_J054613.2-000602':
        ID = "V1647 Ori"
    if ID =='JCMTPP_J053529.8-045944':
        ID = "HOPS 383"
    if ID =='JCMTPP_J162626.8-242431':
        ID = "VLA 1623-243"
    if ID =='JCMTPP_J182949.8+011520':
        ID = "SMM 1"
    if ID =='JCMTPP_J182948.2+011644':
        ID = "SH 2-68 N"
    if ID =='JCMTPP_J182951.2+011638':
        ID = "EC 53"
    if ID =='JCMTPP_J182952.0+011550':
        ID = "SMM 10"
    if ID =='JCMTPP_J183004.0-020306':
        ID = "CARMA 7"
    if ID =='JCMTPP_J182937.8-015103':
        ID = "IRAS 18270-0153"
    return(ID)

        
        
        
import importlib
class Call_Brights:
    
    import Table_Reader as TR
    importlib.reload(TR)
    def __init__(self):
        input_dir = "/mnt/Secdrive/TSmain/TSmain/"
        rgs = ['IC348','NGC1333','NGC2024','NGC2068','OMC23','OPHCORE','SERPM','SERPS']
        rra  =[56.075, 52.225, 85.42083, 86.55417, 83.87917, 246.77083, 277.45417, 277.50833]
        rdec =[32.08306, 31.28111, -1.8975, -0.10139, -5.01056, -24.54361, 1.25556, -2.04667]
        
        self.F_mean_80 = []
        self.fid_err_80 = []
        self.chi2_80 = []
        self.whole_mean = []
        self.fluxes = []
        
        for i in range(0,8):
            data = TR.ReadTable(input_dir+'Data'+'/'+rgs[i]+'_source_data.dat')
            self.fluxes.append(data.pfluxes)
            self.proto_dists = [np.float16(x[4]) for x in self.a]
            self.disk_dists = [np.float16(x[5]) for x in self.a]
            Type = "X"
            if np.min([data.proto_dists[j], data.disk_dists[j]]) < 10:
                if data.proto_dists[j] < data.disk_dists[j]: Type = "P" 
                else: Type = "D"
            indices = [int(x) for x in data.Indices]
            for j in indices:
                Dist = np.sqrt((rra[i]-data.RAs[j])**2+(rdec[i]-data.DECs[j])**2)*60*60
                if data.mean_pfluxes[j] > 0.014*10 and Dist < 1080:
                    num = len(data.pfluxes[j])
                    num10 = int(np.round(num/10))
                    data_sort = sorted(data.pfluxes[j])
                    data_minmax = np.sum(data_sort[:num10] + data_sort[-num10:])
                    self.F_mean_80.append((data.mean_pfluxes[j]*num-data_minmax)/(num-2*num10))
                    self.fid_err_80.append(np.sqrt(0.014**2+(0.02*F_mean_80[-1])**2))
                    self.chi2_80.append(np.sum(abs(np.array(data_sort[num10:-num10])-F_mean_80[-1])**2/fid_err_80[-1]**2)/(num-2*num10))
                    self.whole_mean.append(data.mean_pfluxes[j])
            

        return(0)
        

class All_Table:
    # Read the table generated after running the main cell of "All_table_plotter.ipynb"
    def __init__(self,inputfile,hexcoord=0):
        self.inputfile = inputfile
        self.a = np.loadtxt(self.inputfile, dtype='str')
        self.Indice = [int(x[0]) for x in self.a]
        self.Region = [x[1] for x in self.a]
        self.IDs = [x[2] for x in self.a]
        if not hexcoord:
            self.RA= [float(x[3]) for x in self.a]
            self.DEC = [float(x[4]) for x in self.a]
        else:
            self.RA= [x[3] for x in self.a]
            self.DEC = [x[4] for x in self.a]
        self.Dist = [float(x[5]) for x in self.a]
        self.Type = [x[6] for x in self.a]
        self.F_mean = [float(x[7]) for x in self.a]
        self.chi2_null = [float(x[8]) for x in self.a]
        self.P_null = [float(x[9]) for x in self.a]
        self.Slope = [float(x[10]) for x in self.a]
        self.dslope = [float(x[11]) for x in self.a]

        self.chi2_Lin = [float(x[13]) for x in self.a]
        #self.P_Lin = [float(x[14]) for x in self.a]
        self.FAP_Lin = [float(x[14]) for x in self.a]
        self.Mean_Sin = [float(x[15]) for x in self.a]
        self.Period = [float(x[16]) for x in self.a]
        self.Amp = [float(x[17]) for x in self.a]
        self.Iphase = [float(x[18]) for x in self.a]
        self.chi2_Sin = [float(x[19]) for x in self.a]
        #self.P_Sin = [float(x[21]) for x in self.a]
        #self.Power_Sin = [float(x[22]) for x in self.a]
        self.FAP_Sin = [float(x[20]) for x in self.a]
        self.modFAP_Sin = [float(x[21]) for x in self.a]
        self.F_amp = [float(x[22]) for x in self.a]
        self.Freq_min = [float(x[23]) for x in self.a]
        self.Freq_max = [float(x[24]) for x in self.a]
        self.dfreq = [float(x[25]) for x in self.a]
        self.F_max = [float(x[26]) for x in self.a]
        self.F_min = [float(x[27]) for x in self.a]
        self.STCH = [float(x[28]) for x in self.a]
        #self.SD_MU = [float(x[29]) for x in self.a]
        
    def AddSpitzer(self):
        # Add Spitzer source data
        # Input: Spitzer source matched table
        self.Tbol = [np.float32(x[56]) for x in self.a]
        self.Lbol = [np.float32(x[57]) for x in self.a]
        #self.Type = [np.float
        
    def AddHOPS(self):
        # Add Herchel Orion Protostars source data
        # Input: HOPS source matched table
        self.Tbol = [float(x[34]) for x in self.a]
        self.Lbol = [float(x[33]) for x in self.a]

class HOPStab:
    
    def __init__(self,inputfile):
        self.inputfile = inputfile
        self.a = np.loadtxt(self.inputfile, dtype='str')
        self.IDs= [x[0] for x in self.a]
        self.RA= [float(x[1]) for x in self.a]
        self.DEC = [float(x[2]) for x in self.a]
        self.Lbol = [float(x[4]) for x in self.a]
        self.Tbol = [float(x[5]) for x in self.a]

        
def Cutto80(data,mean=1,test=0):
    # Cut the top 10 and bottom 10% from the data
    import numpy as np
    num = len(data)
    num10 = int(np.round(num/10))
    data_sort = sorted(data)
    data_minmax = np.sum(data_sort[:num10] + data_sort[-num10:])
    #mean_temp=np.mean(data)-data_minmax
    F_mean_temp = (np.sum(data)-data_minmax)/(num-2*num10)
    #print(mean_temp)
    if mean:
        F_mean_80 =mean+F_mean_temp
    if test:
        print([np.mean(data),F_mean_80,data_minmax])
    fid_err_80 = np.sqrt(0.014**2+(0.02*F_mean_80)**2)
    #chi2_80 = np.sum(abs(np.array(data_sort[num10:-num10])-F_mean_80)**2/fid_err_80**2)/(num-2*num10)
    chi2_80 = np.sum(abs(np.array(data_sort[num10:-num10])-F_mean_temp)**2/fid_err_80**2)/(num-2*num10)
    sigma = abs(data/fid_err_80)/np.max([np.sqrt(chi2_80),1])
    wheremax = np.where(sigma == np.max(sigma))
    #print([sigma, wheremax[0][0]])
    a = data/fid_err_80/np.max([np.sqrt(chi2_80),1])
    Sig_max =a[wheremax[0][0]]

    return(F_mean_80, fid_err_80, chi2_80, Sig_max, wheremax[0][0])
    #num = len(data.pfluxes[j])
    #num10 = int(np.round(num/10))
    #data_sort = sorted(data.pfluxes[j])
    #data_minmax = np.sum(data_sort[:num10] + data_sort[-num10:])
    #F_mean_80.append((data.mean_pfluxes[j]*num-data_minmax)/(num-2*num10))
    #fid_err_80.append(np.sqrt(0.014**2+(0.02*F_mean_80[-1])**2))
    #chi2_80.append(np.sum(abs(np.array(data_sort[num10:-num10])-F_mean_80[-1])**2/fid_err_80[-1]**2)/(num-2*num10))
    #whole_mean.append(data.mean_pfluxes[j])
            
def Columnread(data):
    output = []
    for i in range(0,len(data[0])): 
        output.append([x[i] for x in data])
    return(output)


