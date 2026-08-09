"""
Created on Thu Aug  6 22:17:27 2026

@author: raj1294
"""
#Reduce XMM (EPIC-pn and EPIC-MOS) data and extract time averaged spectra
import numpy as np
from astropy.io import fits
import glob, os
import argparse

ks = 1000
loc = os.getcwd() + "/"

def arguments():

    def bool_to_str(value):
        
        if(value=='True' or value=='False'):
            value = value
        if(value!='True' and value!='False'):
            raise Exception("Enter either True or False")
            
        return value
    
    def str_to_str(value):
            
        value = str(value)
            
        return value
    
    def dm_method(value):
        
        value = str(value)
            
        if(value!='HEA' and value!='SAS'):
            raise Exception("Enter either HEA or SAS")
        
        return value
    
    def egrid_str(valarr):
        
        valarr = valarr.split(",")
        method = valarr[0]
        
        if(method=="minmax"):
            
            emin = float(valarr[1])
            emax = float(valarr[2])
            Nenergy = int(valarr[3])
            egrid = np.logspace(np.log10(emin),np.log10(emax),Nenergy)
            
        if(method=='array'):
            
            elist = valarr[1:]
            for kl in range(len(elist)):
                elist[kl] = float(elist[kl])
            egrid = np.array(elist)
            
        return egrid

    #Arguments to code
    parser = argparse.ArgumentParser(\
    description='Generate images in a specified energy-band and use '+\
    'this to detect sources by cross-correlating with Mexican-Hat functions '+\
    'and generate time-averaged energy spectra'+\
    ' (N.B. Requires and installation of a CIAO environment using conda')

    parser.add_argument('-srad','--srcrad',default=0.0083333,\
    help='Source extraction radius [in degrees]',\
    required=True,type=float)
        
    parser.add_argument('-brad','--bkgrad',default=0.0083333,\
    help='Background extraction radius [in degrees]',\
    required=True,type=float)
    
    parser.add_argument('-refemin','--reference_energy_min',default=0.3,\
    help='Reference-band minimum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-refemax','--reference_energy_max',default=12.0,\
    help='Reference-band maximum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-srcdet','--sigthresh',default=1e-16,\
    help='Source detection threshold',\
    required=True,type=float)

    args = vars(parser.parse_args())
    
    return args

# Pipeline to construct lags spectra from XMM catalogue
# Source detection parameters
args = arguments()

#Reference-band minimum energy [in keV]
refemin = args['reference_energy_min']
#Reference-band maximum energy [in keV]
refemax = args['reference_energy_max']

#Source detection threshold
sigthresh = args['sigthresh']
#Source extraction region size
srcradius = args['srcrad']
#Background extraction region size
srcradiusbkg = args['bkgrad']

#Separation between catalogue ULX position and source detection 
septhresh = 0.04 
dsep = 1e-5 #Adaptive separation step
pimin = refemin*1000 #Min channel of image
pimax = refemax*1000 #Max channel of image

#Keys
dirkey = "0*/proc/source_list*epn*.fits" 
obsidkey = "0*"

ctref = 0
stringdet = []
stringdet.append("mkdir PSFs/")
stringdet.append("mkdir lags/")
for ObsId in sorted(glob.glob(obsidkey)):    
    ctr = 1
    
    for evim in sorted(glob.glob(ObsId + "/proc/*EPN*Imaging*.ds")):
                        
        hdu = fits.open(evim)
        telapse = hdu[1].header['ONTIME']
        
        unfiltfile = "epn_obs" + str(ctr) + ".fits"
        filtfile = "epnclean_obs" + str(ctr) + ".fits"
        
        commdet1 = "cp " + evim + " " + loc + ObsId + "/proc/" + unfiltfile
        commdet2 = "cd " + ObsId + "/proc/"
                    
        #Remove background flares
        bkgrate = "rateepn" + str(ctr) + ".fits"
        gtifile = "gtiepn" + str(ctr) + ".fits"
        
        commdet3 = "evselect table=" + unfiltfile +\
        " energycolumn=PI " +\
        "expression='#XMMEA_EP && (PATTERN==0)" +\
        " && (PI>10000 && PI<12000)'" +\
        " withrateset=yes rateset=" + bkgrate + " timebinsize=100" +\
        " maketimecolumn=yes makeratecolumn=yes"
        
        commdet4 = "RT=$(bkgoptrate tssettabname=" + bkgrate +\
        " | sed -n 3p)"    
        
        commdet5 = 'tabgtigen table=' + bkgrate +\
                ' expression="RATE <= $RT" ' + "gtiset=" + gtifile
                  
        commdet6 = "evselect table=" + unfiltfile + " withfilteredset=Y" +\
                " filteredset=" + filtfile +\
                " destruct=Y keepfilteroutput=T" +\
                ' expression="#XMMEA_EP && (FLAG==0)' +\
                ' && (PATTERN==0) &&' +\
                ' gti(' + gtifile + ',TIME) && (PI>=' + str(pimin) +\
                ' && PI<=' + str(pimax) + ')" '
        
        stringdet.append(commdet1)
        stringdet.append(commdet2)
        stringdet.append(commdet3)
        stringdet.append(commdet4)
        stringdet.append(commdet5)
        stringdet.append(commdet6)
        
        #Create image
        imfile = "epnimage" + str(ctr) + ".fits"
        commdet7 = "evselect table=" + filtfile +\
        " withimageset=yes imageset=" +\
        imfile + " xcolumn=X ycolumn=Y" +\
        " imagebinning=imageSize ximagesize=600 yimagesize=600"
        imfilenew = "epnimage_" + ObsId + "_" + str(ctr) + ".fits"
        
        commdet7b = "cp " + imfile + " ../../PSFs/" + imfilenew
        stringdet.append(commdet7)
        stringdet.append(commdet7b)
        
        #Run source detection algorithm
        srcfile = "source_list_epn" + str(ctr) + ".fits"
        commdet8 = "source activate ciao-4.17"
        commdet9 = "wavdetect " + imfile +\
        " source_list_epn.fits source_cell.fits image.fits" +\
        " background.fits expfile=none psffile=none " +\
        'scales="1 2 4 8 16" sigthresh=' + str(sigthresh) +\
        ' regfile=. clobber=yes'
        commdet10 = "mv source_list_epn.fits " + srcfile
        commdet11 = "echo $PWD"
        commdet12 = "conda deactivate"
        commdet13 = "cd ../../"
        commdet14 = ""
        
        stringdet.append(commdet8)
        stringdet.append(commdet9)
        stringdet.append(commdet10)
        stringdet.append(commdet11)
        stringdet.append(commdet12)
        stringdet.append(commdet13)
        stringdet.append(commdet14)
        
        ctr += 1

np.savetxt("srcdetect_epn.sh",stringdet,fmt='%s',delimiter='   ')
os.system("chmod u+x *.sh")
os.system("./srcdetect_epn.sh")
