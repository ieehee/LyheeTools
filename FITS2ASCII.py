from astropy.io import fits
from astropy.wcs import WCS
import os
import pathlib
import numpy as np
import matplotlib.pyplot  as plt


def Get_fitinfo(filename, path="./"):
    image = fits.open(path+filename)
    data = image[0].data
    header = image[0].header
    RA_size = header["NAXIS1"]
    Dec_size = header["NAXIS2"]
    RA_Refpix = header["CRPIX1"]
    Dec_Refpix = header["CRPIX2"]
    RA_Refdeg = header["CRVAL1"]
    Dec_Refdeg = header["CRVAL2"]
    RA_Deldeg = header["CDELT1"]
    Dec_Deldeg = header["CDELT2"]
    
    RAs = np.linspace(RA_Refdeg-RA_Refpix*RA_Deldeg, RA_Refdeg+RA_Refpix*RA_Deldeg, RA_size)
    Decs = np.linspace(Dec_Refdeg-Dec_Refpix*Dec_Deldeg, Dec_Refdeg+Dec_Refpix*Dec_Deldeg, Dec_size)
    return(data, RAs, Decs)


def FITS2ASCII(RAs, Decs, data, filename="./OUTPUT.ascii"):
    with open(filename, 'w') as w:
        w.write("#RA Dec Flux")
        for i, RA in enumerate(RAs):
            for j, Dec in enumerate(Decs):
                if len(data.shape) == 2:
                    output = '\n'+str(RA)+' '+str(Dec)+' '+str(data[j,i])
                else:
                    output = '\n'+str(RA)+' '+str(Dec)+' '+str(data[0,j,i])
                w.write(output)
    print('File saved at "'+filename+'"')
    return(0)

#Example below%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls_dir = os.listdir()
fitlist = [x for x in ls_dir if x.find(".fit")+1] # +1 because find method return -1 if there is nothing.

for file in fitlist:
    data, RAs, Decs = Get_fitinfo(file)
    FITS2ASCII(RAs, Decs, data, filename=file[:-4]+".ascii")
