#Produce background subtracted light-curves
import numpy as np
from astropy.io import fits
import glob, os
from astropy.coordinates import SkyCoord
import argparse

ks = 1000

#Set current location directory
loc = os.getcwd() + "/"

# ObsID(s)
obsid = "0*"

#Arguments to code
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
    
    parser = argparse.ArgumentParser(\
    description='Generate light-curves in a'+\
                ' reference energy band and comparison energy-bands to extract '+\
                'covariance spectra or for pulsation and QPO searches')
        
    parser.add_argument('-srcname','--sourcename',default='',\
    help='Target Name',required=True,type=str_to_str)

    parser.add_argument('-srcrad','--sourceradius',default='',\
    help='Target Name',required=True,type=str_to_str)

    parser.add_argument('-dtbincov','--binning_time_cov',default=100,\
    help='Binning time of LC for covariance spectra [in seconds]',\
    required=True,type=float)

    parser.add_argument('-dtbinqpo','--binning_time_qpo',default=100,\
    help='Binning time of LC for QPO searches [in seconds]',\
    required=True,type=float)

    parser.add_argument('-texp','--threshold_exp_time',default=5,\
    help='Only use observation if elapsed time > tthresh [in ks]',\
    required=True,type=float)

    parser.add_argument('-sthresh','--source_det_thresh',default=1e-16,\
    help='Source detection threshold',\
    required=True,type=float)
        
    parser.add_argument('-psearch','--pulse_search',default='False',\
    help='Extract light-curves with DT=DTMIN?',\
    required=True,type=bool_to_str)

    parser.add_argument('-bsub','--bkg_sub',default='False',\
    help='Use epiclccorr to subtract background?',\
    required=True,type=bool_to_str)

    parser.add_argument('-aflag','--add_flag',default='False',\
    help='Add Flag to extract events?',\
    required=True,type=bool_to_str)

    parser.add_argument('-refemin','--reference_energy_min',default=0.3,\
    help='Reference-band minimum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-refemax','--reference_energy_max',default=12.0,\
    help='Reference-band maximum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-egrid','--energy_grid',default='False',\
    help='Specify comparison-band energy grid [in keV]',\
        required=True,type=egrid_str)
    
    parser.add_argument('-rmflares','--remove_bkg_flares',default='False',\
    help='Remove background flares?',required=True,type=bool_to_str)

    ags = vars(parser.parse_args())
    
    return ags

args = arguments()

#Search for pulsations or QPOs?
psearch = args['pulse_search']
#Subtract background?
bkgsubepiclc = args['bkg_sub']
addflag = args['add_flag']
#Separation between target source location and ULX location (in degrees) 
srcradius = args['sourceradius'] #source extraction region size (in degrees)
#LC bin time (for cov spec)
bintimecov = args['binning_time_cov'] 
#Set to frame time (for pulsations or QPOs)
bintimepulse = args['binning_time_qpo'] 
#Threshold exposure time
tthresh = args['threshold_exp_time']*ks
Emin = args['reference_energy_min'] #Ref EMIN
Emax = args['reference_energy_max'] #Ref EMAX
#Source detection
srcname = args['sourcename']
csrc = SkyCoord.from_name(srcname)
ra_src = csrc.ra.degree
dec_src = csrc.dec.degree

#Threshold exposure time
rmflares = args['remove_bkg_flares']

#Energy grid
energies = args['energy_grid']
pimin = 1000*np.array(energies[0:-1])
pimax = 1000*np.array(energies[1:])
pirefmin = Emin*1000
pirefmax = Emax*1000

#Fixed parameters
septhresh = 0.04 #Threshold separation
dsep = 1e-5 #Adaptive separation step

stringcov,stringpulse = [[],[]]

for ObsId in sorted(glob.glob(obsid+"*")):
    
    ObsId = ObsId.split("/")[0]
                                                                                    
    ctr = 1
    for evim in sorted(glob.glob(loc + ObsId + "/proc/*EPN*Imaging*.ds")):
                                                                                                                                
        unfiltfile = "epn_obs" + str(ctr) + ".fits"
        unfiltfilebkg = "epn_bkg" + str(ctr) + ".fits"
        srcfile = loc + ObsId + "/proc/source_list_epn" + str(ctr) + ".fits" 
        visnum = srcfile.split("/")[-1].split("_")[-1].split(".fits")[0].\
                 split("epn")[1]
                
        commv1 = "cp " + evim + " " + loc + ObsId + "/proc/" +\
        unfiltfile
        commv2 = "cp " + evim + " " + loc + ObsId + "/proc/" +\
        unfiltfilebkg
                
        hdu = fits.open(srcfile)
        tbdata = hdu[1].data
        RA = tbdata['RA']
        DEC = tbdata['DEC']        
        cts = tbdata['NET_COUNTS']
        sigma = tbdata['SRC_SIGNIFICANCE']
                
        hdutime = fits.open(evim)
        telapse = hdutime[1].header['ONTIME']
                                                                
        ctr += 1
                
        culx = SkyCoord(ra_src,dec_src,unit="deg",frame='icrs')                            
        csrc = SkyCoord(RA,DEC,unit="deg",frame="icrs")
                    
        ra = csrc.ra.degree
        dec = csrc.dec.degree   
        separation = culx.separation(csrc)
        sepdeg = separation.deg
                                                                    
        ranew = ra[sepdeg<=septhresh]
        decnew = dec[sepdeg<=septhresh]
        ctsnew = cts[sepdeg<=septhresh]
        sigmanew = sigma[sepdeg<=septhresh]
        sepdegnew = sepdeg[sepdeg<=septhresh]
        
        bkgrate = "rateepn_obs" + str(srcfile[-6]) + ".fits"
        gtifile = "gtiepn_obs" + str(srcfile[-6]) + ".fits"
        unfiltfileclean = "epn_obs" + str(srcfile[-6]) + ".fits"
        unfiltfilecleanbkg = "epn_bkg" + str(srcfile[-6]) + ".fits"
        bkgrate_new = "rateepn_obs" + ObsId + str(srcfile[-6]) +\
                      ".fits"
                                                                                                                                
        if(len(ranew)>1):
                            
            while(True):
                
                ranew = ranew[sepdegnew<=septhresh]
                decnew = decnew[sepdegnew<=septhresh]
                ctsnew = ctsnew[sepdegnew<=septhresh]
                sigmanew = sigmanew[sepdegnew<=septhresh]
                sepdegnew = sepdegnew[sepdegnew<=septhresh]
                septhresh -= dsep
                
                if(septhresh<0):
                    
                    dsep *= 0.1
                    break
                
                if(len(ranew)==1):
                    break
        
        if(len(ranew)==1  and telapse>tthresh):
                                                        
            radet = ranew[0]
            decdet = decnew[0]
            ctsdet = ctsnew[0]
            sigdet= sigmanew[0]
            sepdet = sepdeg[0]
            
            stringcov.append(commv1)
            stringcov.append(commv2)
            stringpulse.append(commv1)
            stringpulse.append(commv2)
            
            commmv3 = "cd " + ObsId + "/proc"
            stringcov.append(commmv3)
            stringpulse.append(commmv3)
                                            
            #Update CCF directory
            commccf = "export SAS_CCF=" + loc + ObsId +\
                      "/proc/ccf.cif"
            stringcov.append(commccf)
            stringpulse.append(commccf)
            
            #Update ODF directory
            for sasfile in glob.glob(loc + ObsId + "/proc/*.SAS"):
                commodf = "export SAS_ODF=" + sasfile
                stringcov.append(commodf)
                stringpulse.append(commodf)
            
            #Barycenter event file
            unfiltfile = "epn_obs" + str(srcfile[-6]) + ".fits"
            commbary = "barycen table=" + unfiltfile + ":EVENTS " +\
            "timecolumn=TIME withsrccoordinates=yes srcra=" +\
                      str(radet) + " srcdec=" + str(decdet) +\
                      " ephemeris=DE405"
            stringcov.append(commbary)
            stringcov.append("")

##############################################################################
            
            #Remove background flares
            if(rmflares=="True"):
                                              
                commbg1 = "evselect table=" + unfiltfile +\
                         " expression='#XMMEA_EP && (PATTERN==0) &&" +\
                         " (PI>10000 && PI<12000)' withrateset=yes" +\
                         " rateset=" + str(bkgrate) + " timebinsize=" +\
                         str(bintimecov) +\
                         " maketimecolumn=yes makeratecolumn=yes"
                
                commbg2 = "cp " + str(bkgrate) + " ../../lags/" +\
                            bkgrate_new
                commbg3 = "ratethresh=$(bkgoptrate tssettabname=" +\
                             bkgrate + " | sed -n 3p)"
                commbg4 = "tabgtigen table=" + bkgrate +\
                           ' expression="RATE <= ' +\
                           str('$ratethresh"') + " gtiset=" + str(gtifile)
                
                commbg5 = "evselect table=" + unfiltfile +\
                           " withfilteredset=Y" +\
                           " filteredset=" + unfiltfileclean +\
                " destruct=Y keepfilteroutput=T expression='#XMMEA_EP &&" +\
                " (FLAG==0) && (PATTERN<=4) && gti" +\
                "(" + str(gtifile) + ",TIME) && (PI>=300 && PI<=12000)'"
            
                stringcov.append(commbg1)
                stringcov.append(commbg2)
                stringcov.append(commbg3)
                stringcov.append(commbg4)
                stringcov.append(commbg5)
                stringpulse.append(commbg1)
                stringpulse.append(commbg2)
                stringpulse.append(commbg3)
                stringpulse.append(commbg4)
                stringpulse.append(commbg5)

#############################################################################
            
            #Filtered event file (full-band)
            filtevrefst = "epn_src_obs_" + ObsId +\
                          "_ref.fits"
            newfiltevrefst = "epn_net_obs_" + ObsId + "_" +\
                             str(srcfile[-6]) + "_ref.fits"
                             
            commfullband = "evselect table=" + unfiltfileclean +\
            " expression='#XMMEA_EP && (FLAG==0)" +\
            " && (PATTERN<=4) && PI in [300:12000]" +\
            " && (RA,DEC) in CIRCLE(" + str(radet) + "," +\
            str(decdet) + "," + str(srcradius) +\
            ")' withfilteredset=Y filteredset=" +\
            filtevrefst    
            
            stringcov.append(commfullband)
            stringcov.append("")
            
#############################################################################
            
            #Comparison band LCs
            for k in range(len(energies)-1):
                                    
                PIMIN = str(int(pimin[k]))
                PIMAX = str(int(pimax[k]))
                
                pirefmin1 = str(int(pirefmin))
                pirefmax1 = PIMIN
                pirefmin2 = PIMAX
                pirefmax2 = str(int(pirefmax))
                                                                                                                            
                filtlcref = "epn_src_obs" + ObsId + "_" +\
                            str(srcfile[-6]) +\
                            "_en" + str(k+1) + "_ref.lc"
                filtlcbkg = "epn_bkg_obs" + ObsId + "_" +\
                            str(srcfile[-6]) +\
                            "_en" + str(k+1) + "_ref.lc" 
                newlcref = "epn_net_obs" + ObsId + "_" +\
                            str(srcfile[-6]) +\
                           "_en" + str(k+1) + "_ref.lc"
                filtbkgref = "epn_bkg_obs" + ObsId + "_" +\
                             str(srcfile[-6]) +\
                             "_en" + str(k+1) + "_ref.fits"                    
                filtevref = "epn_src_obs" + ObsId + "_" +\
                            str(srcfile[-6]) +\
                            "_en" + str(k+1) + "_ref.fits"
                                        
                #Subtract channel of interest
                if(k!=0 and k!=len(energies)-2):
                    
                    strcommadd = "PI in [" + str(pirefmin1) + ":" +\
                                 str(pirefmax2) + "]" +\
                                 " && " +\
                                 "(!PI in [" + str(int(pimin[k])) +\
                                 ":" + str(int(pimax[k])) + "])"
                    
                if(k==0):
                    
                    strcommadd = "(PI in [" + str(pirefmin2) +\
                    ":" + str(pirefmax2) + "])"

                if(k==len(energies)-2):
                    
                    strcommadd = "(PI in [" + str(pirefmin1) +\
                    ":" + str(pirefmax1) + "])"
                
                #Filtered reference-band light-curve
                commrefbandlc = "evselect table=" + unfiltfileclean +\
                         " expression='#XMMEA_EP && (FLAG==0)" +\
                         " && (PATTERN<=4) && " +\
                         str(strcommadd) + " && " +\
                         "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                         str(decdet) + "," + str(srcradius) + ")'" +\
                         " rateset=Y rateset=" +\
                         filtlcref + " maketimecolumn=Y timebinsize=" +\
                         str(bintimecov) + " makeratecolumn=Y"     
                stringcov.append(commrefbandlc)

                #Filtered reference-band event file
                newfiltev_ref = "epn_net_obs_" + ObsId + "_" +\
                              str(srcfile[-6]) +\
                              "_en" + str(k+1) + "_ref.fits"
                              
                commrefband = "evselect table=" + unfiltfileclean +\
                              " expression='#XMMEA_EP && (FLAG==0)" +\
                              " && (PATTERN<=4) && " +\
                              str(strcommadd) + " && " +\
                           "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                           str(decdet) + "," + str(srcradius) +\
                           ")' withfilteredset=Y filteredset=" +\
                           filtevref    
                
                #Move filtered reference-band event file
                stringcov.append(commrefband)                    

###############################################################################

                #Background extraction region
                with open(ObsId + "/proc/bkg1.reg") as fi:
                    
                    for line in fi:
                        line = line.split()
                        if(len(line[0])>5):
                            if(line[0][0:6]=='circle'):
                        
                                raback = float(line[0].split(",")[0][7:])
                                decback = float(line[0].split(",")[1])
                                                                
                                #Filtered background event file
                                commbackrefev =\
                                "evselect table=" + unfiltfileclean +\
                                " expression='#XMMEA_EP && (FLAG==0)" +\
                                " && (PATTERN<=4) && " + str(strcommadd) +\
                                " && (RA,DEC) in CIRCLE(" + str(raback) +\
                                "," + str(decback) + "," +\
                                str(srcradius) +\
                                ")' withfilteredset=Y filteredset=" +\
                                filtbkgref 
                                stringcov.append(commbackrefev)
                                                                         
                                #Filtered background light-curve
                                commbacklc = "evselect table=" +\
                                unfiltfileclean +\
                                " expression='#XMMEA_EP && (FLAG==0)" +\
                                " && (PATTERN<=4) && " +\
                                str(strcommadd) +\
                                " && (RA,DEC) in CIRCLE(" + str(raback) +\
                                "," +\
                                str(decback) + "," + str(srcradius) +\
                                ")'" +\
                                " rateset=Y rateset=" + filtlcbkg +\
                                " maketimecolumn=Y timebinsize=" +\
                                str(bintimecov) + " makeratecolumn=Y"
                                stringcov.append(commbacklc)
                                                                         
                                #Background subtract (reference-band)
                                if(bkgsubepiclc=="True"):
                                    
                                    commbkgsub = "epiclccorr srctslist=" +\
                                    filtlcref + " eventlist=" +\
                                    unfiltfileclean + " outset=" +\
                                    newlcref +\
                                    " bkgtslist=" + filtlcbkg +\
                                    " withbkgset=yes" +\
                                    " applyabsolutecorrections=yes"
                                    stringcov.append(commbkgsub)
                                    
###############################################################################

                #Broadband LC
                filtlcrefst = "epn_src_obs_" + ObsId +\
                              "_ref.lc"
                filtlcbkgrefst = "epn_bkg_obs_" + ObsId +\
                              "_ref.lc"
                newfiltlcrefst = "epn_net_obs_" + ObsId +\
                                 "_ref.lc"
                
                if(k==0 and psearch=="True"):
                    
                                            
                    #Barycenter source event file
                    unfiltfile = "epn_obs" + str(srcfile[-6]) + ".fits"
                    commbarysrc = "barycen table=" + unfiltfile + ":EVENTS " +\
                    "timecolumn=TIME withsrccoordinates=yes srcra=" +\
                              str(radet) + " srcdec=" + str(decdet) +\
                              " ephemeris=DE405"
                    
                    #Barycenter background event file
                    unfiltfilebkg = "epn_bkg" + str(srcfile[-6]) + ".fits"
                    commbarybkg = "barycen table=" + unfiltfilebkg +\
                    ":EVENTS " +\
                    "timecolumn=TIME withsrccoordinates=yes srcra=" +\
                    str(raback) + " srcdec=" + str(decback) +\
                    " ephemeris=DE405"
                    
                    stringpulse.append(commbarysrc)
                    stringpulse.append(commbarybkg)
                                            
                    #Source and background light-curves
                    if(addflag=="False"):
                        
                        commsrclc = "evselect table=" + unfiltfileclean +\
                                   " expression='#XMMEA_EP" +\
                                   " && (PATTERN<=4)" +\
                                   " && PI in [300:12000] && " +\
                                   "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                                   str(decdet) + "," + str(srcradius) + ")'" +\
                                   " rateset=Y rateset=" +\
                                   filtlcrefst +\
                                   " maketimecolumn=Y timebinsize=" +\
                                   str(bintimepulse) + " makeratecolumn=Y"     
                        
                        commbkglc = "evselect table=" + unfiltfilecleanbkg +\
                                   " expression='#XMMEA_EP" +\
                                   " && (PATTERN<=4)" +\
                                   " && PI in [300:12000] && " +\
                                   "(RA,DEC) in CIRCLE(" +\
                                   str(raback) + "," +\
                                   str(decback) + "," + str(srcradius) +\
                                   ")'" +\
                                   " rateset=Y rateset=" +\
                                   filtlcbkgrefst +\
                                   " maketimecolumn=Y timebinsize=" +\
                                   str(bintimepulse) + " makeratecolumn=Y" 
                    
                        #Filtered event file
                        filtev = "epn_src_obs" + ObsId + "_" +\
                                 str(srcfile[-6]) + ".fits"
                        newfiltev = "epn_net_obs_" + str(srcfile[-6]) +\
                                    "_" + ObsId + "_ref.fits"                         
                        commfullev = "evselect table=" + unfiltfile +\
                                     " withfilteredset=Y filteredset=" +\
                                     filtev +\
                        " expression='#XMMEA_EP" +\
                        " && (PATTERN<=4) && (PI in [300:12000]) && " +\
                        "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                        str(decdet) + "," + str(srcradius) + ")'"
                                                
                    
                    if(addflag=="True"):
                        
                        #Source and background light-curves
                        commsrclc = "evselect table=" + unfiltfileclean +\
                                   " expression='#XMMEA_EP && (FLAG==0)" +\
                                   " && (PATTERN<=4)" +\
                                   " && PI in [300:12000] && " +\
                                   "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                                   str(decdet) + "," + str(srcradius) + ")'" +\
                                   " rateset=Y rateset=" +\
                                   filtlcrefst +\
                                   " maketimecolumn=Y timebinsize=" +\
                                   str(bintimepulse) + " makeratecolumn=Y"     
                        
                        commbkglc = "evselect table=" + unfiltfilecleanbkg +\
                                   " expression='#XMMEA_EP && (FLAG==0)" +\
                                   " && (PATTERN<=4)" +\
                                   " && PI in [300:12000] && " +\
                                   "(RA,DEC) in CIRCLE(" +\
                                   str(raback) + "," +\
                                   str(decback) + "," + str(srcradius) +\
                                   ")'" +\
                                   " rateset=Y rateset=" +\
                                   filtlcbkgrefst +\
                                   " maketimecolumn=Y timebinsize=" +\
                                   str(bintimepulse) + " makeratecolumn=Y"  
                
                        #Filtered event file
                        filtev = "epn_src_obs" + ObsId + "_" +\
                                 str(srcfile[-6]) + ".fits"
                        newfiltev = "epn_net_obs_" + str(srcfile[-6]) +\
                                    "_" + ObsId + "_ref.fits"                         
                        commfullev = "evselect table=" + unfiltfile +\
                                     " withfilteredset=Y filteredset=" +\
                                     filtev +\
                        " expression='#XMMEA_EP && (FLAG==0)" +\
                        " && (PATTERN<=4) && (PI in [300:12000]) && " +\
                        "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                        str(decdet) + "," + str(srcradius) + ")'"
                        
                    stringpulse.append(commsrclc)
                    stringpulse.append(commbkglc)
                    stringpulse.append(commfullev)
                    
                    if(bkgsubepiclc=="False"):
                                                    
                        stringpulse.append("mv " + filtlcrefst +\
                        " ../../lags/" + filtlcrefst)
                        stringpulse.append("mv " + filtlcbkgrefst +\
                        " ../../lags/" + filtlcbkgrefst)
                        stringpulse.append("mv " + filtev +\
                        " ../../lags/" + newfiltev)

                    if(bkgsubepiclc=="True"):
                                                
                        stringpulse.append("mv " + filtlcrefst +\
                        " ../../lags/" + newfiltlcrefst)
                        stringpulse.append("mv " + filtlcbkgrefst +\
                        " ../../lags/" + filtlcbkgrefst)
                                            
                #Filtered comparison-band lightcurve
                srclc_comp = "epn_src_obs" + ObsId + "_" +\
                             str(srcfile[-6]) +\
                             "_en" + str(k+1) + "_comp.lc"
                commfiltlccomp = "evselect table=" + unfiltfileclean +\
                               " expression='#XMMEA_EP && (FLAG==0)" +\
                          " && (PATTERN<=4) && (PI in [" + PIMIN +\
                          ":" + PIMAX + "]) && " +\
                          "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                          str(decdet) + "," + str(srcradius) + ")'" +\
                          " rateset=Y rateset=" +\
                          srclc_comp +\
                          " maketimecolumn=Y timebinsize=" +\
                          str(bintimecov) + " makeratecolumn=Y"
                stringcov.append(commfiltlccomp)

                #Filtered comparison-band event file
                filtev_comp = "epn_src_obs" + ObsId + "_" +\
                              str(srcfile[-6]) +\
                              "_en" + str(k+1) + "_comp.fits"
                newfiltev_comp = "epn_net_obs" + ObsId + "_" +\
                              str(srcfile[-6]) +\
                              "_en" + str(k+1) + "_comp.fits"
                commfiltevcomp = "evselect table=" + unfiltfileclean +\
                          " withfilteredset=Y filteredset=" +\
                          filtev_comp +\
                          " expression='#XMMEA_EP && (FLAG==0)" +\
                          " && (PATTERN<=4) && (PI in [" + PIMIN +\
                          ":" + PIMAX + "]) && " +\
                          "(RA,DEC) in CIRCLE(" + str(radet) + "," +\
                          str(decdet) + "," + str(srcradius) + ")'"
                stringcov.append(commfiltevcomp)
                
                #Filtered background comparison-band event file
                filtbkgcomp = "epn_bkg_obs" + ObsId + "_" +\
                              str(srcfile[-6]) +\
                              "_en" + str(k+1) + "_comp.fits"     
                              
                commbkgevcomp =\
                "evselect table=" + unfiltfileclean +\
                " expression='#XMMEA_EP && (FLAG==0)" +\
                " && (PATTERN<=4) && " + "(PI in [" + PIMIN +\
                ":" + PIMAX + "]) && " +\
                " (RA,DEC) in CIRCLE(" + str(raback) +\
                "," + str(decback) + "," +\
                str(srcradius) +\
                ")' withfilteredset=Y filteredset=" +\
                filtbkgcomp
                stringcov.append(commbkgevcomp)
                    
                #Filtered background comparison-band lightcurve
                filtlcbkg_comp = "epn_bkg_obs" + ObsId + "_" +\
                                 str(srcfile[-6]) + "_en" +\
                                 str(k+1) + "_comp.lc"               
                commbkgcomplc = "evselect table=" + unfiltfileclean +\
                         " expression='#XMMEA_EP && (FLAG==0)" +\
                         " && (PATTERN<=4) && (PI in [" + PIMIN +\
                         ":" + PIMAX + "]) && " +\
                         "(RA,DEC) in CIRCLE(" + str(raback) + "," +\
                         str(decback) + "," + str(srcradius) + ")'" +\
                         " rateset=Y rateset=" + filtlcbkg_comp +\
                         " maketimecolumn=Y timebinsize=" +\
                         str(bintimecov) + " makeratecolumn=Y"
                stringcov.append(commbkgcomplc)
                
                #Subtract background (comparison-band)
                newlc_comp = "epn_net_obs" + ObsId + "_" +\
                             str(srcfile[-6]) + "_en" +\
                             str(k+1) + "_comp.lc"
                newlcbkgref = "epn_bkg_obs" + ObsId + "_" +\
                            str(srcfile[-6]) + "_en" +\
                            str(k+1) + "_ref.lc"
                newlcbkgcomp = "epn_bkg_obs" + ObsId + "_" +\
                                str(srcfile[-6]) + "_en" +\
                                str(k+1) + "_comp.lc"
                
                #Reference-band and comparison-band LCs  
                if(bkgsubepiclc=="False"):

                    stringcov.append("mv " + filtlcref +\
                                  " ../../lags/" +\
                                  newlcref)
                    stringcov.append("mv " + srclc_comp +\
                                  " ../../lags/" +\
                                  newlc_comp)
                    stringcov.append("mv " + filtlcbkg +\
                                  " ../../lags/" +\
                                  newlcbkgref)
                    stringcov.append("mv " + filtlcbkg_comp +\
                                  " ../../lags/" +\
                                  newlcbkgcomp)
                
                #Background reference-band and comparison-band LCs  
                if(bkgsubepiclc=="True"):
                                        
                    commbkgepiclc = "epiclccorr srctslist=" + srclc_comp +\
                              " eventlist=" +\
                              unfiltfileclean + " outset=" +\
                              newlc_comp +\
                              " bkgtslist=" + filtlcbkg_comp +\
                              " withbkgset=yes" +\
                              " applyabsolutecorrections=yes"
                    stringcov.append(commbkgepiclc)
                                                                
                    stringcov.append("mv " + newlcref +\
                                     " ../../lags/" +\
                                     newlcref)
                    stringcov.append("mv " + newlc_comp +\
                                     " ../../lags/" +\
                                     newlc_comp)
                    stringcov.append("mv " + newlcbkgref +\
                                     " ../../lags/" +\
                                     newlcbkgref)
                    stringcov.append("mv " + newlcbkgcomp +\
                                     " ../../lags/" +\
                                     newlcbkgcomp)

                    stringcov.append("")

            
            stringcov.append("cd ../../")
            stringcov.append("")
            stringcov.append("")
            
            stringpulse.append("cd ../../")
            stringpulse.append("")
                    
np.savetxt("filter_lc_cov_epn.sh",stringcov,fmt='%s',delimiter='   ')
os.system("chmod u+x filter_lc_cov_epn.sh")
os.system("./filter_lc_cov_epn.sh")

# np.savetxt("filter_lc_pulse.sh",stringpulse,fmt='%s',delimiter='   ')
# os.system("chmod u+x filter_lc_pulse.sh")
# os.system("./filter_lc_pulse.sh")

