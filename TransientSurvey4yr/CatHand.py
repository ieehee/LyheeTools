####Catalogue Handler#######

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import importlib

class LoadDunham08:
    def __init__(self):
        
        file = "/mnt/Secdrive/TSmain/TSmain/VeLLOs/Dunham08.cat"
        T,Group, indice, RAh, RAm, RAs, DECd, DECm, DECs, LIRs=np.loadtxt(file,unpack=True,dtype='str')
        self.Indice = indice
        self.RAh = [float(x) for x in RAh]
        self.RAm = [float(x) for x in RAm]
        self.RAs = [float(x) for x in RAs]
        self.DECd = [float(x) for x in DECd]
        self.DECm = [float(x) for x in DECm]
        self.DECs = [float(x) for x in DECs]
        self.LIRs = [float(x) for x in LIRs]

        
#class LoadKim19:
#    def __init__(self):
#        file = "/mnt/Secdrive/TSmain/TSmain/VeLLOs/Kim19.cat"
#        #ID, RA, DEC, Region, Dist, L_int, M_env, Class =np.loadtxt(file,unpack=True,dtype='str',delimiter='\t')
#        self.test =np.loadtxt(file,unpack=True,dtype='str',delimiter='\t')
#        #aaa= aaa
#        #self.Indice = indice
#        self.RAh = [float(np.split(x,' ')[0]) for x in RA]
#        self.RAm = [float(np.split(x,' ')[1]) for x in RA]
#        self.RAs = [float(np.split(x,' ')[2]) for x in RA]
#        self.DECd = [float(np.split(x,' ')[0]) for x in DEC]
#        self.DECm = [float(np.split(x,' ')[1]) for x in DEC]
#        self.DECs = [float(np.split(x,' ')[2]) for x in DEC]
#        self.L_int = [x for x in L_int]

        
        
class LoadKim16:
    def __init__(self):
        file = "/mnt/Secdrive/TSmain/TSmain/VeLLOs/Kim16.cat"
        indices, IDs, Regions, RAs, DECs, Lints, e_Lints, Tbols, Lbols, T_dusts, types =np.loadtxt(file,unpack=True,dtype='str')
        self.Indice  = [float(x) for x in indices]
        self.IDs = [x for x in IDs]
        self.Regions = [x for x in Regions]
        self.RA = [float(x) for x in RAs]
        self.DEC = [float(x) for x in DECs]
        self.Lints = [float(x) for x in Lints]
        self.e_Lints = [float(x) for x in e_Lints]
        self.Tbols = [float(x) for x in Tbols]
        self.Lbols = [float(x) for x in Lbols]
        self.T_dusts = [float(x) for x in T_dusts]
        self.types = types

class LoadWISE:
    def __init__(self):
        file = '/mnt/Secdrive/TSmain/TSmain/Data/WISE_source.tab'
        x, ind,dist, RAs, DECs, W1_mean, W1_stdev, W1_emean, W2_mean, W2_stdev, W2_emean =np.loadtxt(file,unpack=True, dtype='str')
        #x =np.loadtxt(file,unpack=True, dtype='str')
        self.index = [float(x) for x in ind]
        self.RA = [float(x) for x in RAs]
        self.DEC = [float(x) for x in DECs]

class LoadWISEvar:
    def __init__(self):
        file = '/mnt/Secdrive/TSmain/TSmain/Data/WISE_var.tab'
        x, ind,dist, RAs, DECs, W1_mean, W1_stdev, W1_emean, W2_mean, W2_stdev, W2_emean,  vartype , classes=np.loadtxt(file,unpack=True, dtype='str')
        #x =np.loadtxt(file,unpack=True, dtype='str')
        self.index = [float(x) for x in ind]
        self.RA = [float(x) for x in RAs]
        self.DEC = [float(x) for x in DECs]
        self.vartype = [x for x in vartype]
        self.classes = [x for x in classes]
        
class LoadWISEproto:
    def __init__(self):
        file = '/mnt/Secdrive/TSmain/TSmain/Data/WISE_protos.tab'
        X=np.loadtxt(file, dtype='str')
        #x =np.loadtxt(file,unpack=True, dtype='str')
        self.index = [float(x[1]) for x in X]
        self.RA = [float(x[3]) for x in X]
        self.DEC = [float(x[4]) for x in X]
        #self.vartype = [x[12] for x in X]
        self.classes = [x[11] for x in X]
        

class ALMAcodes:
    def __init__(self):
        file = '/mnt/Secdrive/TSmain/TSmain/ALMA_Arxiv/ALMA_Codes.tsv'
        X= np.loadtxt(file, dtype='str')
        self.Pcode = [x[0] for x in X]
        self.ALMAID = [x[1] for x in X]
        self.RA = [float(x[1]) for x in X]
        self.DEC = [float(x[2]) for x in X]
        self.Band = [float(x[3]) for x in X]
        self.Release = [x[4] for x in X]
        self.Pubs = [float(x[5]) for x in X]
        self.Ang_Res = [float(x[6]) for x in X]
        self.MaxRC = [float(x[7]) for x in X]
        self.FOV = [float(x[8]) for x in X]
        
class CompareCats:
    def __init__(self, coord1, coord2,unit=0):
        #self.Indice1,self.IDs1,self.RAs1,self.DECs1 = np.loadtxt(file1,dtype='str')
        #self.Indice2,self.IDs2,self.RAs2,self.DECs2 = np.loadtxt(file1,dtype='str')
        
        c1 = SkyCoord(coord1,frame='icrs')
        c2 = SkyCoord(coord2,frame='icrs')
        if unit:
            c1 = SkyCoord(coord1,frame='icrs',unit=unit)
            c2 = SkyCoord(coord2,frame='icrs',unit=unit)
        print(file1 + ' Loaded')
        print(file2 + ' Loaded')
        
        
        