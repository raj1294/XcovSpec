#Reduce XMM (EPIC-pn and EPIC-MOS) data and extract time averaged spectra
import numpy as np
import glob, os
import argparse
from astropy.coordinates import SkyCoord
from astroquery.heasarc import Heasarc

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
    
    def obsid_str(valarr):
                
        valarr = valarr.split(",")
        obslist = valarr[:]
        obsids = np.array(obslist)
            
        return obsids

    #Arguments to code
    parser = argparse.ArgumentParser(\
    description='Reduce XMM-Newton data using epproc and emproc, '+\
    'generate images in a specified energy-band and use '+\
    'this to detect sources using Mexican-Hat correlation method '+\
    'and generate time-averaged energy spectra')

    parser.add_argument('-srcname','--sourcename',\
    help='Target Name',required=True,type=str_to_str)

    parser.add_argument('-obsids','--observation_ids',default='',\
    help='List of Observation Identifiers',\
    required=True,type=obsid_str)
        
    parser.add_argument('-refemin','--reference_energy_min',default=0.3,\
    help='Reference-band minimum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-refemax','--reference_energy_max',default=12.0,\
    help='Reference-band maximum energy [keV]',\
    required=True,type=float)

    args = vars(parser.parse_args())
    
    return args

# Pipeline to construct lags spectra from XMM catalogue
# Source detection parameters
args = arguments()
#Source name
srcname = args['sourcename']
#Reference-band minimum energy [in keV]
refemin = args['reference_energy_min']
#Reference-band maximum energy [in keV]
refemax = args['reference_energy_max']
pimin = refemin*1000 #Min channel of image
pimax = refemax*1000 #Max channel of image
#Observation identifiers
obsidsel = args['observation_ids']
#Separation between catalogue ULX position and source detection 
septhresh = 0.04 
dsep = 1e-5 #Adaptive separation step

#Source name
pos = SkyCoord.from_name(srcname)
tab = Heasarc.query_region(pos, catalog='xmmmaster')
obsid = tab['obsid'].value

#Download and Reduce
stringspec,stringreduce,stringdet,stringdownload = [[],[],[],[]]

#Obs IDs
ObsId = []
for jsel in range(len(obsidsel)):    
    for ksel in range(len(obsid)):
        if(obsidsel[jsel]==obsid[ksel]):     
            ObsId.append(obsidsel[jsel])
ObsId = np.array(ObsId)

for kobs in range(len(ObsId)):
                        
    comm1 = "mkdir " + ObsId[kobs]
    comm2 = "mkdir " + ObsId[kobs]  + "/odf/"
    
    stringdownload.append(comm1)
    stringdownload.append(comm2)
            
    outfile = str(ObsId[kobs]) + ".tar"
    
    url_download =\
    "'https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno="+\
    str(ObsId[kobs])+"&level=ODF'"   
    
    comm3 = "wget -O " + outfile + " " + url_download
    stringdownload.append(comm3)
        
    comm4 = "mv *.tar " + ObsId[kobs] + "/odf/"        
    comm5 = "cd " + ObsId[kobs] + "/odf/"
    comm6 = "tar -xvf *.tar"
    comm7 = "tar -xvf *.TAR"
    comm8 = "rm -f *.tar *.TAR"
    comm9 = "cd ../"
    comm10 = "mkdir proc/"
    comm11 = "cd ../"
                        
    stringdownload.append(comm4)
    stringdownload.append(comm5)
    stringdownload.append(comm6)
    stringdownload.append(comm7)
    stringdownload.append(comm8)
    stringdownload.append(comm9)
    stringdownload.append(comm10)
    stringdownload.append(comm11)
    stringdownload.append("")

np.savetxt("download.sh",stringdownload,fmt='%s',delimiter='  ')
os.system("chmod u+x download.sh")
os.system("./download.sh")

#Requires the installation of SAS
for mobs in range(len(ObsId)):
    
    comm0 = "cd " + ObsId[mobs] + "/"
    comm1 = "z=$(pwd)"
    comm2 = "cd proc/"
    comm3 = "export SAS_CCF=$z" + "/odf"
    comm4 = "export SAS_ODF=$z" + "/odf"
    comm5 = "cifbuild"
    comm6 = "z=$(pwd)"
    comm7 = "export SAS_CCF=$z" + "/ccf.cif"
    comm8 = "odfingest"
    comm9 = "q=(*.SAS)"
    comm10 = "export SAS_ODF=$z" + "/$q"
    comm11 = "epproc"
    comm12 = "emproc"
    comm13 = "cd ../../"
    comm14 = ""
            
    stringreduce.append(comm0)
    stringreduce.append(comm1)
    stringreduce.append(comm2)
    stringreduce.append(comm3)
    stringreduce.append(comm4)
    stringreduce.append(comm5)
    stringreduce.append(comm6)
    stringreduce.append(comm7)
    stringreduce.append(comm8)
    stringreduce.append(comm9)
    stringreduce.append(comm10)
    stringreduce.append(comm11)
    stringreduce.append(comm12)
    stringreduce.append(comm13)
    stringreduce.append(comm14)
    
np.savetxt("reduce_xmm.sh",stringreduce,fmt='%s',delimiter='   ')
os.system("chmod u+x *.sh")
os.system("./reduce_xmm.sh")