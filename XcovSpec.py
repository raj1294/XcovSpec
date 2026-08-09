import numpy as np
import argparse
import os

ks = 1000

stringcomms = []

def arguments():
    
    def str_to_str(value):
            
        value = str(value)
            
        return value

    def bool_to_psd(valarr):
            
        #Raise Exception
        if(valarr!='abs' and valarr!='frac' and valarr!='leahy'):
            raise Exception("Enter either abs, frac or leahy")
        return valarr

    def float_to_str(valarr):
        
        valarrsplit = valarr.split(",")
        gbin,gbinstingray = valarrsplit[0],valarrsplit[1]
        return gbin, gbinstingray

    def bool_to_str(value):
        
        if(value=='True' or value=='False'):
            value = value
        if(value!='True' and value!='False'):
            raise Exception("Enter either True or False")
            
        return value
          
    def bool_to_str_seg(valarr):
        
        valarrsplit = valarr.split(",")
        
        if(valarrsplit[0]=='True'):
            valarrsplit[1] = float(valarrsplit[1])
            valarrsplit[2] = float(valarrsplit[2])

        #If not using segment LC, use the entire LC
        if(valarrsplit[0]=='False'):
            valarrsplit[1] = 0
            valarrsplit[2] = 1e90
        
        #Raise Exception
        if(valarrsplit[0]!='True' and valarrsplit[0]!='False'):
            raise Exception(\
            "Enter either True or False for first argument")
            
        return valarrsplit[0],valarrsplit[1],valarrsplit[2]

    def bool_to_str_gaps(valarr):
        
        valarrsplit = valarr.split(",")
        
        #Raise Exception
        if(valarrsplit[0]!='True' and valarrsplit[0]!='False'):
            raise Exception(\
            "Enter either True or False for first argument")
        
        if(valarrsplit[0]=='True'):
            
            if(valarrsplit[1]!='B' and valarrsplit[1]!='T' and\
               valarrsplit[1]!='S'):
                raise Exception(\
                "Enter either B, T or S for second argument")
            
        return valarrsplit[0],valarrsplit[1]

    def bool_to_str_mcmc(valarr):
        
        valarrsplit = valarr.split(",")
        
        #Raise Exception
        if(valarrsplit[0]!='True' and valarrsplit[0]!='False'):
            raise Exception(\
            "Enter either True or False for first argument")
        
        if(valarrsplit[0]=='False'):
            valarrsplit[1] = 0
            
        return valarrsplit[0],valarrsplit[1]

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
            
        return egrid, valarr

    def obsid_str(valarr):
        
        valarr = valarr.split(",")
        obslist = valarr[:]
        obsids = np.array(obslist)
            
        return obsids

    #Arguments to code
    parser = argparse.ArgumentParser(\
    description='Generate covariance spectra, lag-energy spectra '+\
                'plot power-spectral densities and plot light-curves')

    parser.add_argument('-srcname','--sourcename',\
    help='Target Name',required=True,type=str_to_str)

    parser.add_argument('-plc','--plotlc',\
    help='Plot LC? [Enter either True or False]',\
    required=True,type=bool_to_str,default=True)
                        
    parser.add_argument('-plags','--plotlags',\
    help='Plot Lag-energy spectrum? [Enter either True or False]',\
    required=True,type=bool_to_str,default=True)
        
    parser.add_argument('-ppsd','--plotpsd',\
    help='Plot PSD? [Enter either True or False]',\
    required=True,type=bool_to_str,default=True)
        
    parser.add_argument('-split','--splitscheme',\
    help="Segment LC? [Enter either True or False]",\
    required=True,type=bool_to_str)

    parser.add_argument('-statpower','--statspsd',\
    help='Does the LC follow a Poisson distribution? ' +\
    '[Enter either True or False]',required=True,type=bool_to_str,\
    default=True)

    parser.add_argument('-normpower','--normpsd',\
    help='Enter normalisation of PSD [abs or frac or leahy]',\
    required=True,type=bool_to_psd,default=True)    

    parser.add_argument('-refemin','--reference_energy_min',default=0.3,\
    help='Reference-band minimum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-refemax','--reference_energy_max',default=12.0,\
    help='Reference-band maximum energy [keV]',\
    required=True,type=float)

    parser.add_argument('-flgaps','--fillgaps',\
    help='Fill LC gaps? [Enter 2 values separated by comma of type: '+\
    'Boolean(Enter either True or False)'+\
    ' String(Interpolation scheme: Enter B (bootstrapping)'+\
    ' ,T (timmer-koenig) or S (window deconvolution))',\
    required=True,type=bool_to_str_gaps)
                        
    parser.add_argument('-seglc','--segmentlc',\
    help='Segment LC? [Enter 3 values separated by comma of type: ' +\
    'True/False (Boolean), Start of LC in ks (float), ' +\
    'End of LC in ks (float)]',\
    required=True,type=bool_to_str_seg)

    parser.add_argument('-fmin','--freqmin',default=1e-4,\
    help='Minimum Fourier Frequency for Covariance '+\
    '[Enter a floating point value in Hz]',\
    required=True,type=float)
        
    parser.add_argument('-fmax','--freqmax',default=5e-4,\
    help='Maximum Fourier Frequency for Covariance '+\
    '[Enter a floating point value in Hz]',\
    required=True,type=float)
        
    parser.add_argument('-gbin','--geombin',\
    help='Geometric binning factor [Enter two floating point values]',\
    required=True,type=float_to_str)
        
    parser.add_argument('-gscale','--groupscale',default=1,\
    help=\
    'Scaling factor to group Covariance spectrum [Enter an integer value]',\
    required=True,type=int)

    parser.add_argument('-rmcmc','--runmcmc',\
    help='Run MCMC? [Enter 2 values separated by comma of type: '+\
    'Boolean(Enter either True or False)'+\
    ', Int(Enter number of MCMC simulations)',\
    required=True,type=bool_to_str_mcmc,default=True)

    parser.add_argument('-psdmods','--powspecmod',\
    help="Model power spectral density? [Enter either True or False]",\
    required=True,type=bool_to_str)
    
    parser.add_argument('-egrid','--energy_grid',default='False',\
    help='Specify comparison-band energy grid [in keV]',\
    required=True,type=egrid_str)

    parser.add_argument('-gencov','--covspec',\
    help="Generate Covariance Spectrum? [Enter either True or False]",\
    required=True,type=bool_to_str)

    parser.add_argument('-mincts','--minimum_cts',default=30,\
    help='Minimum counts in source spectrum',\
    required=True,type=int)

    parser.add_argument('-minbcts','--minimum_cts_bkg',default=30,\
    help='Minimum counts in background spectrum',\
    required=True,type=int)

    parser.add_argument('-srad','--srcrad',default='',\
    help='Target Name',required=True,type=str_to_str)

    parser.add_argument('-brad','--bkgrad',default=0.0083333,\
    help='Background extraction radius [in degrees]',\
    required=True,type=float)

    parser.add_argument('-rmflares','--remove_bkg_flares',default='False',\
    help='Remove background flares?',required=True,type=bool_to_str)

    parser.add_argument('-srcdet','--sigthresh',default=1e-16,\
    help='Source detection threshold',\
    required=True,type=float)

    parser.add_argument('-dtbinbkg','--bkg_bin_time',default=100,\
    help='Background binning time [in seconds]',\
    required=True,type=int)

    parser.add_argument('-dtbincov','--binning_time_cov',default=100,\
    help='Binning time of LC for covariance spectra [in seconds]',\
    required=True,type=float)

    parser.add_argument('-aflag','--add_flag',default='False',\
    help='Add Flag to extract events?',\
    required=True,type=bool_to_str)

    parser.add_argument('-obsids','--observation_ids',default='',\
    help='List of Observation Identifiers',\
    required=True,type=obsid_str)

    parser.add_argument('-texp','--threshold_exp_time',default=1,\
    help='Only use observation if elapsed time > tthresh [in ks]',\
    required=False,type=float)

    parser.add_argument('-dtbinqpo','--binning_time_qpo',default=100,\
    help='Binning time of LC for QPO searches [in seconds]',\
    required=True,type=float)

    parser.add_argument('-psearch','--pulse_search',default='False',\
    help='Extract light-curves with DT=DTMIN?',\
    required=True,type=bool_to_str)

    parser.add_argument('-bsub','--bkg_sub',default='False',\
    help='Use epiclccorr to subtract background?',\
    required=True,type=bool_to_str)

    ags = vars(parser.parse_args())
    
    return ags

args = arguments()

# PSD parameters 
gbinning = args['geombin']
gbinningargs = str(gbinning[0]) + "," + str(gbinning[1])

freqmin = args['freqmin']
freqmax = args['freqmax']
plotpsd = args['plotpsd']
statpow = args['statspsd']
splitscheme = args['splitscheme']
normpsd = args['normpsd']
psdmods = args['powspecmod']
gencov = args['covspec']

#Reference band energy (min and max)
Emin = args['reference_energy_min']
Emax = args['reference_energy_max']

#LC parameters
plotlc = args['plotlc']
fillgaps,fillmethod = args['fillgaps']
segmentlc,tmin,tmax = args['segmentlc']

if(fillmethod=="B"):
    bootstrap = "True"
    timmerkoenig = "False"
    stdwin = "False"
if(fillmethod=="T"):
    timmerkoenig = "True"
    bootstrap = "False"
    stdwin = "False"
if(fillmethod=="S"):
    stdwin = "True"
    bootstrap = "False"
    timmerkoenig = "False"

#Generate covspec
gencov = args['covspec']

#Source radius
srcrad = args['srcrad']
bkgrad = args['bkgrad']

#MCMC simulations
rmcmc = args['runmcmc']
rmcmcargs = str(rmcmc[0]) + ',' + str(rmcmc[1])

#Plot lags
plotlags = args['plotlags']

#Group covariance spectrum
groupscale = args['groupscale']

#Energy grid
energies = args['energy_grid'][0]
stringargs = args['energy_grid'][1][:]
energyargs = ''
for jen in range(len(stringargs)):
    if(jen!=len(stringargs)-1):
        energyargs += stringargs[jen] + ','
    if(jen==len(stringargs)-1):
        energyargs += stringargs[jen]
        
#Obs ID
obsidsel = args['observation_ids'][0]

#Minimum source counts
mincts = args['minimum_cts']

#Minimum background counts
minctsbkg = args['minimum_cts_bkg']

#Source name
srcname = args['sourcename']
srcname = "'" + srcname + "'"

#Detection threshold 
sigthreshold = args['sigthresh']

#Remove flares
removeflares = args['remove_bkg_flares']

#Background bin time
bkgflarebintime = args['bkg_bin_time']

#Exposure time
tthresh = args['threshold_exp_time']

#Bin time (CV spec and QPQs)
bintimecov = args['binning_time_cov'] 
bintimeqpo = args['binning_time_qpo'] 

#Pulse search
psearch = args['pulse_search']

#Subtract background
bkgsubepiclc = args['bkg_sub']

#Add flag
addflag = args['add_flag']

tmin *= ks
tmax *= ks

commruncovxspec1 =\
'python ~/Documents/SEAWIND/code/Reduction/reduce_xmm.py -srcname ' +\
str(srcname) + ' -obsids ' + str(obsidsel) + ' -refemin ' + str(Emin) +\
' -refemax ' + str(Emax) + ' -mincts ' + str(mincts) + ' -minbcts ' +\
str(minctsbkg) + ' -dtbinbkg ' + str(minctsbkg) 

commruncovxspec2 =\
'python ~/Documents/SEAWIND/code/Reduction/srcdet_xmm.py -srad ' +\
str(srcrad) + ' -brad ' + str(bkgrad) + ' -refemin ' +\
str(Emin) + ' -refemax ' + str(Emax) + ' -srcdet ' +\
str(sigthreshold)

commruncovxspec3 =\
'python ~/Documents/SEAWIND/code/Reduction/extractspec_xmm.py ' +\
' -srcname ' + str(srcname) + ' -rmflares ' + str(removeflares) +\
' -dtbinbkg ' + str(bkgflarebintime) + ' -srad ' +\
str(srcrad) + ' -brad ' + str(bkgrad) + ' -mincts ' +\
str(mincts) + ' -minbcts ' + str(minctsbkg) + ' -refemin ' +\
str(Emin) + ' -refemax ' + str(Emax) + ' -texp ' +\
str(tthresh) + ' -obsids ' + str(obsidsel)

commruncovxspec4 =\
'python ~/Documents/SEAWIND/code/Reduction/reduce_xmm_lc.py -srcname ' +\
str(srcname) + ' -dtbincov ' + str(bintimecov) + ' -dtbinqpo ' +\
str(bintimeqpo) + ' -sthresh ' + str(sigthreshold) +\
' -psearch ' + str(psearch) + ' rmflares ' + str(removeflares) +\
' -bsub ' + str(bkgsubepiclc) + ' -aflag ' + str(addflag) +\
' -refemin ' + str(Emin) + ' -refemax ' + str(Emax) +\
' -egrid ' + str(energyargs) + ' -srcrad ' + str(srcrad) + ' -texp ' +\
str(tthresh)

commruncovxspec5 =\
'python ~/Documents/SEAWIND/code/Timing/covariance.py -plc ' +\
str(plotlc) + ' -plags ' + str(plotlags) + ' -ppsd ' + str(plotpsd) +\
' -split ' + str(splitscheme) + ' -statpower ' +\
str(statpow) + ' -normpower ' + str(normpsd) + ' -flgaps ' +\
str(fillgaps) + ' -seglc ' + str(segmentlc) + ' -fmin ' +\
str(freqmin) + ' -fmax ' + str(freqmax) + ' -gbin ' + str(gbinningargs) +\
' -gscale ' + str(groupscale) + ' -rmcmc ' + str(rmcmcargs) +\
' -gencov ' + str(gencov) + ' -psdmods ' + str(psdmods)

stringcomms.append(commruncovxspec1)
stringcomms.append(commruncovxspec2)
stringcomms.append(commruncovxspec3)
stringcomms.append(commruncovxspec4)
stringcomms.append(commruncovxspec5)
np.savetxt("XcovSpec.sh",stringcomms,fmt='%s',delimiter='  ')
os.system("chmod u+x XcovSpec.sh")
# os.system("./XcovSpec.sh")

