#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 22:21:35 2026

@author: raj1294
"""
#Reduce XMM (EPIC-pn and EPIC-MOS) data and extract time averaged spectra
import numpy as np
from astropy.io import fits
import glob, os
import argparse
from astropy.coordinates import SkyCoord

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

    def obsid_str(valarr):
        
        valarr = valarr.split(",")
        obslist = valarr[:]
        obsids = np.array(obslist)
            
        return obsids

    
    #Arguments to code
    parser = argparse.ArgumentParser(\
    description='Generate time-averaged energy spectra of' +\
    ' reduced XMM observations')
        
    parser.add_argument('-srcname','--sourcename',\
    help='Target Name',required=True,type=str_to_str)
        
    parser.add_argument('-rmflares','--remove_bkg_flares',default='False',\
    help='Remove background flares?',required=True,type=bool_to_str)
    
    parser.add_argument('-dtbinbkg','--bkg_bin_time',default=100,\
    help='Background binning time [in seconds]',\
    required=True,type=int)

    parser.add_argument('-srad','--srcrad',default=0.0083333,\
    help='Source extraction radius [in degrees]',\
    required=True,type=float)
        
    parser.add_argument('-brad','--bkgrad',default=0.0083333,\
    help='Background extraction radius [in degrees]',\
    required=True,type=float)
    
    parser.add_argument('-mincts','--minimum_cts',default=30,\
    help='Minimum counts in source spectrum',\
    required=True,type=int)
    
    parser.add_argument('-minbcts','--minimum_cts_bkg',default=30,\
    help='Minimum counts in background spectrum',\
    required=True,type=int)

    parser.add_argument('-refemin','--reference_energy_min',default=0.3,\
    help='Reference-band minimum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-refemax','--reference_energy_max',default=12.0,\
    help='Reference-band maximum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-obsids','--observation_ids',default='',\
    help='List of Observation Identifiers',\
    required=True,type=obsid_str)

    parser.add_argument('-texp','--threshold_exp_time',default=1,\
    help='Only use observation if elapsed time > tthresh [in ks]',\
    required=False,type=float)

    args = vars(parser.parse_args())
    
    return args

args = arguments()
srcname = args['sourcename'] #Source name

#Source detection
csrc = SkyCoord.from_name(srcname)
ra_src = csrc.ra.degree
dec_src = csrc.dec.degree

#Reference-band minimum and maximum energies [in keV]
refemin = args['reference_energy_min']
refemax = args['reference_energy_max']
pimin = refemin*1000 #Min channel of image
pimax = refemax*1000 #Max channel of image

#Remove bkg flares
removebkgflares = args['remove_bkg_flares'] 
bkgflarebintime = args['bkg_bin_time']

#Source and background extraction region size
srcradius = args['srcrad']
srcradiusbkg = args['bkgrad']

#Minimum number of counts
mincts = args['minimum_cts']
tthresh = args['threshold_exp_time']

#Observation identifiers
obsidsel = args['observation_ids']

#Separation between catalogue ULX position and source detection 
dsep = 1e-5 #Adaptive separation step
septhresh = 0.04 
ctref = 0

if(len(obsidsel)==1):
    dirkey = obsidsel[0] + "/proc/source_list*epn*.fits" 
if(len(obsidsel)!=1):
    dirkey = "0*/proc/source_list*epn*.fits" 
storagedir = "lags/"
stringspec = []
for sourcefile in sorted(glob.glob(loc + dirkey)):
                                                            
    hdu = fits.open(sourcefile)
    tbdata = hdu[1].data
    RA = tbdata['RA']
    DEC = tbdata['DEC']    
    telapse = hdu[1].header['ONTIME']
    ObsId = hdu[1].header['OBS_ID']
    detnum = sourcefile[-6]
                        
    culx = SkyCoord(ra_src,dec_src,unit="deg",frame='icrs')                            
    csrc = SkyCoord(RA,DEC,unit="deg",frame="icrs")
    ranew = csrc.ra.degree
    decnew = csrc.dec.degree   
            
    separation = culx.separation(csrc)
    sepdeg = separation.deg
    ranew = ranew[sepdeg<=septhresh]
    decnew = decnew[sepdeg<=septhresh]
    sepdeg = sepdeg[sepdeg<=septhresh]
        
    if(len(ranew)>1 and telapse>tthresh):
                    
        while(True):
            
            ranew = ranew[sepdeg<=septhresh]
            decnew = decnew[sepdeg<=septhresh]
            sepdeg = sepdeg[sepdeg<=septhresh]
            septhresh -= dsep
                            
            if(len(ranew)==1):
                radet = ranew[0]
                decdet = decnew[0]
                ctref += 1                                            
                break
            
    elif(len(ranew)==1 and telapse>tthresh):
        radet = ranew[0]
        decdet = decnew[0]
        ctref += 1                                            
                    
    commmv = "cd " + ObsId + "/proc"
    stringspec.append(commmv)
    
    #Update CCF directory
    commccf = "export SAS_CCF=" + loc + ObsId +\
              "/proc/ccf.cif"
    stringspec.append(commccf)
    
    #Update ODF directory
    for sasfile in glob.glob(loc + ObsId + "/proc/*.SAS"):
        commodf = "export SAS_ODF=" + sasfile
        stringspec.append(commodf)
    
    stringspec.append("")
    
    unfiltfile = "epn_obs" + str(detnum) + ".fits"
    filtfile = "epn_src" + str(detnum) + ".fits"
    filtfilebkg = "epn_bkg" + str(detnum) + ".fits"
    unfiltfileclean = "epnclean" + str(detnum) + ".fits"            
    unfiltfileclean = unfiltfile       
            
    #Remove background flares
    if(removebkgflares=="yes"):
        
        bkgrate = "rateepn" + str(detnum) + ".fits"
        bkgrate_new = "rateepn" + str(ctref) + ".fits"
        gtifile = "gtiepn" + str(detnum) + ".fits"
        unfiltfileclean = "epnclean" + str(detnum) +\
                          ".fits"            
        unfiltfileclean = unfiltfile           
        commbkgfl1 = "evselect table=" + unfiltfile +\
                     " expression='#XMMEA_EP && (PATTERN==0) &&" +\
                     " (PI>10000 && PI<12000)' withrateset=yes" +\
                     " rateset=" + str(bkgrate) + " timebinsize=" +\
                     str(bkgflarebintime) +\
                   " maketimecolumn=yes makeratecolumn=yes"
        commbkgfl2 = "cp " + str(bkgrate) + " ../../lags/" +\
                    bkgrate_new
        commbkgfl3 = "ratethresh=$(bkgoptrate tssettabname=" +\
                     bkgrate + " | sed -n 3p)"
        commbkgfl4 = "tabgtigen table=" + bkgrate +\
                   ' expression="RATE <= ' +\
                   str('$ratethresh"') + " gtiset=" +\
                   str(gtifile)
        commbkgfl5 = "evselect table=" + unfiltfile +\
        " withfilteredset=Y" + " filteredset=" +\
        unfiltfileclean +\
        " destruct=Y keepfilteroutput=T" +\
        " expression='#XMMEA_EP &&" +\
        " (FLAG==0) && (PATTERN==0) && gti" +\
        "(" + str(gtifile) + ",TIME) && (PI>=300 && PI<=12000)'"
        
        stringspec.append(commbkgfl1)
        stringspec.append(commbkgfl2)
        stringspec.append(commbkgfl3)
        stringspec.append(commbkgfl4)
    
    #Spectrum and response file
    specsrc = "epn_spec" + str(detnum) + ".fits"
    specbkg = "epn_spec" + str(detnum) + "_bkg.fits"
    specrmf = "epn_" + str(ObsId) + "_" + str(detnum) + ".rmf"
    specarf = "epn_" + str(ObsId) + "_" + str(detnum) + ".arf"
    groupspec = "epn_spec" + str(detnum) + "_grp.fits" 
    badpixfile = unfiltfileclean
                  
    #Filtered event file
    commfiltevref = "evselect table=" + unfiltfileclean +\
              " withfilteredset=Y filteredset=" + filtfile +\
              " expression='#XMMEA_EP && (FLAG==0)" +\
              " && (PATTERN==0) && (PI>=300) && (PI<=12000)" +\
              " && " + "(RA,DEC) in CIRCLE(" + str(radet) +\
              "," + str(decdet) + "," + str(srcradius) +\
              ")" + "'"
    commspecbin =  "evselect table=" + filtfile +\
              " withspectrumset=yes " + "spectrumset=" +\
              specsrc + " energycolumn=PI" +\
              " spectralbinsize=5" +\
              " withspecranges=yes specchannelmin=0 " +\
              "specchannelmax=20479"
    commbackscale = "backscale spectrumset=" + specsrc + " " +\
              "badpixlocation=" + badpixfile
    commrmf = "rmfgen spectrumset=" + specsrc + " rmfset=" +\
              specrmf + " extendedsource=no"
    commarf = "arfgen arfset=" + specarf + " spectrumset=" +\
              specsrc + " withrmfset=yes rmfset=" + specrmf +\
              " withbadpixcorr=yes badpixlocation=" +\
              badpixfile + " detmaptype=psf"
                          
    stringspec.append(commfiltevref)
    stringspec.append(commspecbin)
    stringspec.append(commbackscale)
    stringspec.append(commrmf)
    stringspec.append(commarf)
                                
    #Generate background region file
    imagefile = loc + ObsId + "/proc/epnimage" + str(ctref) + ".fits"
    
    print("Please specify a background region file and",\
    "save in proc directory: In order to save it, use the region",\
    "button located in main the panel")
    os.system("ds9 " + imagefile + " -scale log -cmap heat")
        
    for bkgfile in glob.glob(loc + ObsId + "/proc/bkg1.reg"):
                        
        if(len(bkgfile)<=0):
            
            print("Upload background region file")
            
        if(len(bkgfile)>0):
    
            with open(loc + ObsId + "/proc/bkg1.reg") as fi:
                            
                ctbkg = 0
                for line in fi:
                    line = line.split()
                    if(ctbkg==3):
                        raback = float(line[0].split(",")[0][7:])
                        decback = float(line[0].split(",")[1])
                    ctbkg+=1 
                                
                #Filtered background event file
                commbkgev = "evselect table=" + unfiltfileclean +\
                          " withfilteredset=Y filteredset=" +\
                          filtfilebkg +\
                          " expression='#XMMEA_EP && (FLAG==0)" +\
                          " && (PATTERN==0) && " +\
                          "(PI>=300) && (PI<=12000) && " +\
                          "(RA,DEC) in CIRCLE(" + str(raback) + "," +\
                          str(decback) + "," + str(srcradiusbkg) +\
                          ")'"  
                          
                #Pileup file
                pileupfile = "pileup" + str(ObsId) + ".pdf"
                commpileup = "epatplot set=" + str(filtfile) +\
                          " plotfile=" + str(pileupfile) +\
                          " pileupnumberenergyrange='300 12000'" +\
                          " withbackgroundset=yes backgroundset=" +\
                          "" + str(filtfilebkg)
    
                #Background spectrum
                specbkg = "epn_spec" + str(detnum) + "_bkg.fits"
                groupbkg = "epn_spec" + str(detnum) + "_bkg_grp.fits"
                specrmf_bkg = "epn_" + str(ObsId) + "_" + str(detnum) +\
                              "_bkg.rmf"
                specarf_bkg = "epn_" + str(ObsId) + "_" + str(detnum) +\
                              "_bkg.arf"
    
                commbkgspec =  "evselect table=" + filtfilebkg +\
                          " withspectrumset=yes " + "spectrumset=" +\
                          specbkg + " energycolumn=PI " +\
                          "spectralbinsize=5" +\
                          " withspecranges=yes specchannelmin=0 " +\
                          "specchannelmax=20479"
                commback = "backscale spectrumset=" + specbkg + " " +\
                          "badpixlocation=" + badpixfile
                
                commrmf = "rmfgen spectrumset=" + specbkg +\
                            " rmfset=" +\
                          specrmf_bkg + " extendedsource=no"
                commarf = "arfgen arfset=" + specarf_bkg +\
                            " spectrumset=" +\
                            specbkg + " withrmfset=yes rmfset=" +\
                            specrmf_bkg +\
                          " withbadpixcorr=yes badpixlocation=" +\
                          badpixfile + " detmaptype=psf"
                                    
                #Group spectrum                        
                commgrp = "specgroup spectrumset=" + str(specsrc) +\
                          " mincounts=" + str(mincts) +\
                          " oversample=3 rmfset=" + str(specrmf) +\
                          " " + "backgndset=" + str(specbkg) +\
                          " witharfset=yes arfset=" + str(specarf) +\
                          " groupedset=" + str(groupspec)
                commmv = "cd ../../"
                commfill = ""
                
                stringspec.append(commbkgev)
                stringspec.append(commpileup)
                stringspec.append(commbkgspec)
                stringspec.append(commback)
                stringspec.append(commrmf)
                stringspec.append(commarf)
                stringspec.append(commgrp)
                stringspec.append(commmv)
                stringspec.append(commfill)  
                                        
                newdir = ObsId + "/proc/epn*spec*grp.fits"
                spec = "epn_spec" + str(detnum) + "_grp.fits"
                newspec = "epn_spec" + str(detnum) + "_grp_" +\
                ObsId + ".fits"
                bkgspec = "epn_spec" + str(detnum) + "_bkg.fits"
                newbkgspec = "epn_spec" + str(detnum) + "_" +\
                ObsId + "_bkg.fits"     
                newrsp = "epn_" + ObsId + "_" + str(detnum) + ".rmf"
                newarf = "epn_" + ObsId + "_" + str(detnum) + ".arf"
                newspecrmf_bkg = "epn_" + ObsId + "_" + str(detnum) +\
                                 "_bkg.rmf"
                newspecarf_bkg = "epn_" + ObsId + "_" + str(detnum) +\
                                 "_bkg.arf"
                                 
                commmv1 = "cd " + ObsId + "/proc/"
                commmv2 = "mv " + spec + " ../../" + storagedir + newspec
                commmv3 = "mv " + bkgspec + " ../../" +\
                storagedir + newbkgspec
                commmv4 = "mv " + newrsp + " ../../" + storagedir
                commmv5 = "mv " + newarf + " ../../" + storagedir
                commmv6 = "mv " + newspecrmf_bkg + " ../../" + storagedir
                commmv7 = "mv " + newspecarf_bkg + " ../../" + storagedir

                stringspec.append(commmv1)
                stringspec.append(commmv2)
                stringspec.append(commmv3)
                stringspec.append(commmv4)
                stringspec.append(commmv5)
                stringspec.append(commmv6)
                stringspec.append(commmv7)
                stringspec.append("cd ../../")
                stringspec.append(" ")
    

np.savetxt("filter_spec.sh",stringspec,fmt='%s',delimiter='   ')
os.system("chmod u+x filter_spec.sh")
os.system("./filter_spec.sh")

for sourcefile in sorted(glob.glob(loc + dirkey)):

    #Modify exposure time keywords
    newspkey = loc + "lags/epn*spec*grp*" + ObsId + ".fits"
    for newsp in glob.glob(newspkey):
                
        newsp = newsp.split("/")[-1]
        
        hdulistref = fits.open("lags/" + newsp)
        header = hdulistref[2].header
        telapse = header['TELAPSE']
        obsid = newsp.split(".fits")[0].split("grp_")[1]        
        vis = newsp.split("_grp")[0].split("_spec")[1]
        
        newbkg = "epn_spec" + vis + "_" + str(obsid) + "_bkg.fits"
        newrsp = "epn_" + obsid + "_" + vis + ".rmf"
        newarf = "epn_" + obsid + "_" + vis + ".arf"   
        
        commkey1 = "fparkey " + newrsp + " " + storagedir +\
        newsp + "[1] RESPFILE"
        commkey2 = "fparkey " + newarf + " " + storagedir +\
        newsp + "[1] ANCRFILE"
        commkey3 = "fparkey " + newbkg + " " + storagedir +\
        newsp + "[1] BACKFILE"
        
        commkey4 = 'fparkey ' + str(telapse) + " " + storagedir +\
        newsp + '[1] ' + 'EXPOSURE'
        
        commkey5 = 'fparkey ' + str(telapse) + " " + storagedir +\
        newbkg + '[1] ' + 'EXPOSURE'
                
        os.system(commkey1)
        os.system(commkey2)
        os.system(commkey3)
        os.system(commkey4)
        os.system(commkey5)

