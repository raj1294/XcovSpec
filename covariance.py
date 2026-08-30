import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from kapteyn import kmpfit
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.gaussian_process import GaussianProcessRegressor
from stingray import Lightcurve, AveragedPowerspectrum, AveragedCrossspectrum
import warnings, os, glob
from scipy import stats, fft, integrate, special
from astropy.io import fits
from stingray.varenergyspectrum import CovarianceSpectrum
from stingray import EventList
import argparse

warnings.filterwarnings('ignore')

ks = 1000
day = 86400

#Argparse functions
##############################################################################

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
    
    #Arguments to code
    parser = argparse.ArgumentParser(\
    description='Generate covariance spectra, lag-energy spectra '+\
                'plot power-spectral densities and plot light-curves')
    parser.add_argument('-plc','--plotlc',\
    help='Plot LC? [Enter either True or False]',\
    required=True,type=bool_to_str,default=True)

    parser.add_argument('-srcname','--sourcename',default='',\
    help='Target Name',required=True,type=str_to_str)

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

    parser.add_argument('-flgaps','--fillgaps',\
    help='Fill LC gaps? [Enter 2 values separated by comma of type: '+\
    'Boolean(Enter either True or False)'+\
    ' String(Interpolation scheme: Enter B (bootstrapping),T (timmer-koenig) or S (window deconvolution))',\
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
    help='Scaling factor to group Covariance spectrum [Enter an integer value]',\
    required=True,type=int)

    parser.add_argument('-rmcmc','--runmcmc',\
    help='Run MCMC? [Enter 2 values separated by comma of type: '+\
    'Boolean(Enter either True or False), Int(Enter number of MCMC simulations)',\
    required=True,type=bool_to_str_mcmc,default=True)

    parser.add_argument('-gencov','--covspec',\
    help="Generate Covariance Spectrum? [Enter either True or False]",\
    required=True,type=bool_to_str)

    parser.add_argument('-psdmods','--powspecmod',\
    help="Model power spectral density? [Enter either True or False]",\
    required=True,type=bool_to_str)

    ags = vars(parser.parse_args())
    
    return ags

##############################################################################

#Covariance functions

##############################################################################

# Compute Fractional Variability
def Fracvar(rate,error):
    
    mean = np.mean(rate)
    var,varerr = 0,0
    for kfvar in range(len(error)):
        var += (rate[kfvar] - mean)**2
        varerr += error[kfvar]**2
    var/=(len(error)-1)
    varerr/=len(error)
    normexcessvar = (var - varerr)/(mean**2) #Normalised excess variance
    
    Fvar,dFvar = 0,0
    if(normexcessvar>0):
        Fvar = np.sqrt(normexcessvar)
    elif(normexcessvar<0):
        Fvar = 0
        
    errnormexcessvar = np.sqrt(np.sqrt(2.0/len(error))*(varerr/mean**2) +\
    (np.sqrt(varerr/len(error))*(2*Fvar/mean))**2)
    dFvar = np.sqrt(Fvar**2 + errnormexcessvar**2) - Fvar
    
    return Fvar, dFvar

# Window function based on GTIs (Time-Domain)
def rect_window(rate_arr,tS,tE):
    
    twref,rwref = rate_arr        
    
    ywin = np.zeros(len(rwref))
    for lw in range(len(tS)):
        
        gtimin = tS[lw]
        gtimax = tE[lw]
        
        for lw2 in range(len(twref)):
            
            if(twref[lw2]>=gtimin and twref[lw2]<=gtimax):
                
                ywin[lw2] = 1
        
    ywin[0] = 0
    ywin[-1] = 0

    return ywin

# Analytical function to model averaged PSD
def besselmod(pars, xdata):
    
    amplitude,alpha,ampbes,norder = pars
    
    xdata =\
    np.linspace(np.min(xdata),np.max(xdata),len(xdata))/(np.min(xdata))
    ybessel = ampbes*(special.jv(int(norder),xdata))
    ymod = ybessel + amplitude - alpha*xdata 

    return ymod

# Residuals
def resid_besselmod(pars, data):
    
    xdata, ydata, ydataerr = data
    amplitude,alpha,ampbes,norder = pars

    xdata =\
    np.linspace(np.min(xdata),np.max(xdata),len(xdata))/(np.min(xdata))
    resid = (ydata - besselmod(pars,xdata))/ydataerr
    
    return resid 

# Timmer & Koenig (1995) method to generate fake LC
def drawsampbes(freqs,ampfit,expntfit,ampbesfit,norderfit):
    
    ybesselfit = ampbesfit*(special.jv(norderfit,freqs))
    psdomega = ybesselfit + ampfit - expntfit*freqs
    
    r1 = np.random.normal(0.0,scale=np.sqrt(abs(0.5*psdomega)))
    r2 = np.random.normal(0.0,scale=np.sqrt(abs(0.5*psdomega)))
    compnumpos = complex(r1,r2)
    compnumneg = np.conj(compnumpos)
    
    return compnumpos,compnumneg

# Timmer & Koenig (1995) method to generate fake LC
def drawsampgp(psdomegamod):
    
    r1 = np.random.normal(0.0,scale=np.sqrt(abs(0.5*psdomegamod)))
    r2 = np.random.normal(0.0,scale=np.sqrt(abs(0.5*psdomegamod)))
    compnumpos = complex(r1,r2)
    compnumneg = np.conj(compnumpos)
    
    return compnumpos,compnumneg

# Estimation of random noise contribution to lag energy spectrum
def genlc(bin_time_mcmc,tdur,freqmin,freqmax,psdomegamod,murate):
                 
    Nomega = len(psdomegamod)
    complexfft1pos = np.zeros(int(len(psdomegamod))).astype(complex)
    complexfft1neg = np.zeros(int(len(psdomegamod))).astype(complex)

    for irand in range(Nomega):
                        
        compnumber_pos1,compnumber_neg1 = drawsampgp(psdomegamod[irand])
        complexfft1pos[irand] = compnumber_pos1
        complexfft1neg[irand] = compnumber_neg1

        if(irand==Nomega-1):
            complexfft1neg[irand] = np.real(complexfft1neg[irand])
            
    complexfft1neg = np.flip(complexfft1neg)
    complexfft1pos = np.insert(complexfft1pos,0,complex(murate))
    complexfft1 = np.hstack((complexfft1pos,complexfft1neg))
            
    # Artificial LCs generated from PSDs (Timmer & Köenig 1995)
    counts = np.fft.ifft(complexfft1,n=len(complexfft1)) + 0.5*complexfft1[0]
    counts = np.real(counts)
    error = np.sqrt(counts)
    times = bin_time_mcmc*(np.arange(0,len(counts),1))

    # Add a floor
    if(np.min(counts)<0):
        counts -= abs(np.min(counts))
            
    return times, counts, error

#Model PSD with Gaussian Processes
def psdmodgp(tlcpsd,lcpsd,lcerrpsd,reflcpsd,reflcerrpsd,lcbkgpsd,refbkgpsd,\
             Msegpsd,bfactorpsd,Dtpsd,statspsd,rmbtpsd):
           
    # Compute noise level depending on whether the counting statistics
    # are Poissonian or not 
    fnyqpsd = 0.5*(Dtpsd**-1)
    Pnoisepsd,Prefnoisepsd,Msegpsd = 0,0,1
    freqspsd,Pxnpsd,Pynpsd,Cxynpsd,dPxnpsd,dPynpsd,dCxynpsd =\
    [[],[],[],[],[],[],[]]
                
    # Average power spectrum and cross spectrum over M segments
    for kzpsd in range(Msegpsd):
        
        # Split LC into M equal segments
        divpsd = int(len(reflcpsd)/Msegpsd)
        lctemppsd = lcpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
        lctemperrpsd = lcerrpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
        reflctemppsd = reflcpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
        reflctemperrpsd = reflcerrpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
        lcbkgtemppsd = lcbkgpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
        refbkgtemppsd = refbkgpsd[kzpsd*divpsd:(kzpsd+1)*divpsd]
                        
        # Ambient noise level in PSD
        if(statspsd=="True"):
            
            Pnoisepsd += (2*(np.mean(lctemppsd) + np.mean(lcbkgtemppsd))/\
            (np.mean(lctemppsd))**2)
            Prefnoisepsd += (2*(np.mean(reflctemppsd) +\
            np.mean(refbkgtemppsd))/(np.mean(reflctemppsd))**2)
                            
        if(statspsd=="False"):
            
            errsqpsd = 0
            errrefsqpsd = 0
            
            for lpsd in range(len(lctemperrpsd)):
                errsqpsd += lctemperrpsd[lpsd]**2
                errrefsqpsd += reflctemperrpsd[lpsd]**2
                
            errsqpsd/=len(lctemperrpsd)
            errrefsqpsd/=len(reflctemperrpsd)
            Pnoisepsd += errsqpsd/(fnyqpsd*(np.mean(lctemppsd))**2)
            Prefnoisepsd += errrefsqpsd/(fnyqpsd*(np.mean(reflctemppsd))**2) 
        
        if(np.sum(lctemppsd)>0):
            
            lctemppsd = np.array(lctemppsd)
            lctemperrpsd = np.array(lctemperrpsd)
                                                            
            # FFT of comparison-band LC
            Xnpsd = 0.5*(fft.fft(lctemppsd+lctemperrpsd) +\
            fft.fft(lctemppsd-lctemperrpsd))
            Xnerrpsd = 0.5*(fft.fft(lctemppsd+lctemperrpsd)-\
            fft.fft(lctemppsd-lctemperrpsd))
            Xnconjpsd = 0.5*(np.conj(Xnpsd+Xnerrpsd)+\
            np.conj(Xnpsd-Xnerrpsd))
            Xnconjerrpsd = 0.5*(np.conj(Xnpsd+Xnerrpsd)-\
            np.conj(Xnpsd-Xnerrpsd))
            fxnpsd = fft.fftfreq(len(lctemppsd),d=Dtpsd)
                                    
            # FFT of reference-band LC
            Ynpsd = 0.5*(fft.fft(reflctemppsd+reflctemperrpsd)+\
            fft.fft(reflctemppsd-reflctemperrpsd))
            Ynerrpsd = 0.5*(fft.fft(reflctemppsd+reflctemperrpsd)-\
            fft.fft(reflctemppsd-reflctemperrpsd))
            Ynconjpsd = 0.5*(np.conj(Ynpsd+Ynerrpsd)+\
            np.conj(Ynpsd-Ynerrpsd))
            Ynconjerrpsd = 0.5*(np.conj(Ynpsd+Ynerrpsd)-\
            np.conj(Ynpsd-Ynerrpsd))
            fynpsd = fft.fftfreq(len(reflctemppsd),d=Dtpsd)
                                                
            Xnpsd = Xnpsd[fxnpsd>0] 
            Xnerrpsd = Xnerrpsd[fxnpsd>0]
            Xnconjpsd = Xnconjpsd[fxnpsd>0]
            Xnconjerrpsd = Xnconjerrpsd[fxnpsd>0]
            Ynpsd = Ynpsd[fxnpsd>0]
            Ynerrpsd = Ynerrpsd[fxnpsd>0]
            Ynconjpsd = Ynconjpsd[fxnpsd>0]
            Ynconjerrpsd = Ynconjerrpsd[fxnpsd>0]
            fynpsd = fynpsd[fxnpsd>0]
            fxnpsd = fxnpsd[fxnpsd>0]
                                                                                                    
            # Compute PSD and CPSD with rms-squared normalisation for each 
            # segment
            normpsdxpsd =\
            (2.0*Dtpsd)/((len(lctemppsd))*(np.mean(lctemppsd))**2)
            normpsdypsd =\
            (2.0*Dtpsd)/((len(reflctemppsd))*(np.mean(reflctemppsd))**2)
            normcrosspsd =\
            (2.0*Dtpsd)/((len(lctemppsd))*(np.mean(lctemppsd))*\
            (np.mean(reflctemppsd)))
            
            # PSD
            Psdxpsd = normpsdxpsd*Xnconjpsd*Xnpsd
            dPsdxpsd = normpsdxpsd*(Xnconjerrpsd*Xnpsd + Xnconjpsd*Xnerrpsd)
            Psdypsd = normpsdypsd*Ynconjpsd*Ynpsd
            dPsdypsd = normpsdypsd*(Ynconjerrpsd*Ynpsd + Ynconjpsd*Ynerrpsd)
            Crossxypsd = normcrosspsd*Ynconjpsd*Xnpsd
            dCrossxypsd =\
            normcrosspsd*(Ynconjerrpsd*Xnpsd + Ynconjpsd*Xnerrpsd)
                                        
            # Append CPSD and PSDs for each segment to 
            # pass to functions for averaging and binning
            if(len(Crossxypsd)>0 and len(Psdxpsd)>0 and len(Psdypsd)>0):
                                
                freqspsd.append(fxnpsd)
                Pxnpsd.append(Psdxpsd)
                Pynpsd.append(Psdypsd)
                Cxynpsd.append(Crossxypsd)
    
    freqspsd = np.array(freqspsd)
    Pxnpsd = np.array(Pxnpsd)
    dPxnpsd = np.array(dPxnpsd)
    Pynpsd = np.array(Pynpsd)
    dPynpsd = np.array(dPynpsd)
    Cxynpsd = np.array(Cxynpsd)
    dCxynpsd = np.array(dCxynpsd)
    Nfmodrefpsd = 0
                        
    if(len(Cxynpsd)>0):
                        
        # Average PSDs and CPSD over M segements
        freqypsd,Pxavgpsd,dPxavgpsd,Pyavgpsd,dPyavgpsd,\
        Cxyavgpsd,dCxyavgpsd =\
        Pbin(Msegpsd,freqspsd,Pxnpsd,Pynpsd,Cxynpsd)
                
        avgPxpsd = Pxavgpsd
        avgPypsd = Pyavgpsd
        avgCxypsd = Cxyavgpsd
        avgPxerrpsd = dPxavgpsd
        avgPyerrpsd = dPyavgpsd
        avgCxyerrpsd = dCxyavgpsd
        Karrpsd = np.ones(len(Pxavgpsd))
                
        # Implement frequency dependent binning of averaged PSDs and CPSD
        bfactorpsd = 1.1
        if(bfactorpsd>1):
            freqxpsd,avgPxpsd,avgPypsd,avgCxypsd,avgPxerrpsd,\
            avgPyerrpsd,avgCxyerrpsd,Karrpsd =\
            fbin(bfactorpsd,freqypsd,Pxavgpsd,Pyavgpsd,\
                 Cxyavgpsd,dPxavgpsd,dPyavgpsd,dCxyavgpsd)
        
        Karrpsd = np.array(Karrpsd)
        avgPxerrpsd = np.array(avgPxerrpsd)
        avgPyerrpsd = np.array(avgPyerrpsd)
        avgCxyerrpsd = np.array(avgCxyerrpsd)
        avgPxpsd = np.array(avgPxpsd)
        avgPypsd = np.array(avgPypsd)
        avgCxypsd = np.array(avgCxypsd)
        freqxpsd = np.array(freqxpsd)
        freqypsd = freqxpsd
        
        lgfreqxpsd = np.log(freqxpsd)
        lgfreqypsd = np.log(freqypsd)
        lgavgPypsd = np.log(avgPypsd)
        lgavgPyerrpsd = abs(avgPyerrpsd/avgPypsd)
        
        if(len(Cxynpsd)==0):
            
            freqxpsd = np.zeros(Nfmodrefpsd)   
            freqypsd = np.zeros(Nfmodrefpsd)
            avgPxpsd = np.zeros(Nfmodrefpsd)
            avgPxerrpsd = np.zeros(Nfmodrefpsd)
            avgPypsd = np.zeros(Nfmodrefpsd)
            avgPyerrpsd = np.zeros(Nfmodrefpsd)
                    
        # Primary kernel parameters (RBF)
        lscale = 20.0
        sigf = 10
        sign = 4.37e-4
        dim = 1
        lgfreqcompre = lgfreqxpsd.reshape(len(lgfreqxpsd),dim)
        kern = (sigf**2)*RBF(length_scale=lscale) +\
        WhiteKernel(noise_level=sign)
        
        gp = GaussianProcessRegressor(kernel=kern,alpha=1e-10,\
        n_restarts_optimizer=200,normalize_y=True)
        gp.fit(lgfreqcompre,lgavgPypsd)
        
        scorecomp = gp.score(lgfreqcompre,lgavgPypsd)
        paramscomp = gp.kernel_
        
        # Best-fit prediction
        Npsfmod = int(0.5*len(tlcpsd))
        
        lgfreqypsdmod =\
        np.linspace(np.min(lgfreqypsd),np.max(lgfreqypsd),Npsfmod)
        lgfmodcomppsdre = lgfreqypsdmod.reshape(len(lgfreqypsdmod),dim)
        
        lgypsdmod,lgypsdmoderr =\
        gp.predict(lgfmodcomppsdre,return_std=True)
        mulcpsd = np.mean(avgPypsd)
                    
    return lgfreqypsd,lgavgPypsd,lgavgPyerrpsd,lgfreqypsdmod,lgypsdmod,\
           mulcpsd
           
#Estimation of random noise contribution to lag-energy spectrum
def mcmc_det(bin_time_mcmc,tdur,nsegmts,geom_rebin,freqmin,freqmax,\
             A1,ind1,Abes1,norbes1,A2,ind2,Abes2,norbes2,\
             murate1,murate2,plts,sts,method):
             
    #Draw randomly from best-fit PSD and inverse FFT to generate LC
    omegamin = 1./tdur
    omegamax = 0.5*(bin_time_mcmc**-1)
    domega = omegamin
    omega = np.arange(omegamin,omegamax+domega,domega)
    
    complexfft1pos = np.zeros(int(len(omega))).astype(complex)
    complexfft1neg = np.zeros(int(len(omega))).astype(complex)
    complexfft2pos = np.zeros(int(len(omega))).astype(complex)
    complexfft2neg = np.zeros(int(len(omega))).astype(complex)
    
    for irand in range(int(len(omega))):
        
        compnumber_pos1,compnumber_neg1 =\
        drawsampbes(omega[irand],A1,ind1,Abes1,norbes1)
        compnumber_pos2,compnumber_neg2 =\
        drawsampbes(omega[irand],A2,ind2,Abes2,norbes2)
        
        complexfft1pos[irand] = compnumber_pos1
        complexfft1neg[irand] = compnumber_neg1
        complexfft2pos[irand] = compnumber_pos2
        complexfft2neg[irand] = compnumber_neg2

        if(irand==int(len(omega))-1):
            complexfft1neg[irand] = np.real(complexfft1neg[irand])
            complexfft2neg[irand] = np.real(complexfft2neg[irand])
        
    complexfft1neg = np.flip(complexfft1neg)
    complexfft1pos = np.insert(complexfft1pos,0,complex(2*murate1))
    complexfft1 = np.hstack((complexfft1pos,complexfft1neg))
    complexfft2neg = np.flip(complexfft2neg)
    complexfft2pos = np.insert(complexfft2pos,0,complex(2*murate2))
    complexfft2 = np.hstack((complexfft2pos,complexfft2neg))
                
    #Artificial LCs generated from PSDs (Timmer & Köenig 1995)
    counts1 = np.fft.ifft(complexfft1)
    counts1 = np.real(counts1)    
    counts2 = np.fft.ifft(complexfft2)
    counts2 = np.real(counts2)
    
    error1 =\
    np.sqrt(np.random.poisson(np.int64(abs(counts1)*bin_time_mcmc)))/\
    bin_time_mcmc
    error2 =\
    np.sqrt(np.random.poisson(np.int64(abs(counts2)*bin_time_mcmc)))/\
    bin_time_mcmc
    times = np.arange(0,bin_time_mcmc*len(counts1),bin_time_mcmc)
    
    #Ensure positive values for counts while retaining the shape
    ct1_min = abs(np.min(counts1))
    ct1_max = abs(np.min(counts2))
    counts1 = (counts1 + ct1_min)*bin_time_mcmc
    counts2 = (counts2 + ct1_max)*bin_time_mcmc
    error1 = error1*bin_time_mcmc
    error2 = error2*bin_time_mcmc
        
    if(plts=="True"):
        
        plt.errorbar(times,counts1,yerr=error1,fmt='k-')
        plt.errorbar(times,counts2,yerr=error2,fmt='r-')
        plt.show()
    
    #Fake lags
    if(method=="timelags"):
        
        freq_fake,dfreq_fake,phaselag_fake,\
        phaselag_efake,coh_fake,cohe_fake =\
        time_lag_func(counts1,error1,counts2,error2,\
                      np.zeros(len(counts1)),np.zeros(len(counts2)),\
                      np.ones(len(counts1)),nsegmts,geom_rebin,\
                      bin_time_mcmc,sts)
        phaselag_fake = phaselag_fake[freq_fake>=freqmin]
        freq_fake = freq_fake[freq_fake>=freqmin]
        phaselag_fake = phaselag_fake[freq_fake<=freqmax]
        freq_fake = freq_fake[freq_fake<=freqmax]
        lg_fake = phaselag_fake/(2.0*np.pi*freq_fake)
    
    #Fake lags (stingray)
    if(method=="stingray"):
        
        lcref_sim = Lightcurve(times,counts=counts1,err=error1,\
                               dt=bin_time_mcmc)
        lccomp_sim = Lightcurve(times,counts=counts2,err=error2,\
                                dt=bin_time_mcmc)
        evref_sim = EventList.from_lc(lcref_sim)
        evcomp_sim = EventList.from_lc(lccomp_sim)
                
        tsegsize = tdur/nsegmts
        CSAsim = AveragedCrossspectrum.from_events(evref_sim,evcomp_sim,\
                 segment_size=tsegsize,norm="frac",use_common_mean=True,\
                 dt=bin_time_mcmc,silent=True)
        CSAsim = CSAsim.rebin_log(geom_rebin)
        
        freq_fake = CSAsim.freq
        lag_fake, lag_e_fake = CSAsim.time_lag()
        lag_fake = lag_fake[freq_fake>=freqmin]
        freq_fake = freq_fake[freq_fake>=freqmin]
        lag_fake = lag_fake[freq_fake<=freqmax]
        freq_fake = freq_fake[freq_fake<=freqmax]
        lg_fake = np.mean(lag_fake)
    
    return lg_fake

#Remove NANs in lags
def remove_nans_lags(arrnans):
        
    newarrnanslist = []
    for qdnans in range(len(arrnans)):
                                                                              
        isnanarr = np.isnan(arrnans[qdnans])
        for qdnans2 in range(len(isnanarr)):   
            
            if(isnanarr[qdnans2]==True):
                
                newindex = np.arange(0,len(arrnans),1)
                
                for qdnans3 in range(len(newindex)):
                    arrnans[qdnans3][qdnans2] = -1e10
        
        newarrnanslist.append(arrnans[qdnans])

    for pindnan in range(len(newarrnanslist)):    
        newarrnanslist[pindnan] = np.array(newarrnanslist[pindnan])
        newarrnanslist[pindnan] = newarrnanslist[pindnan]\
                                  [newarrnanslist[pindnan]>-1e9]
    newarrnanslist = np.array(newarrnanslist)
    return newarrnanslist

#Remove NANs in LC
def remove_nans_lc(arrnans):
        
    newarrnanslist = []
    for qdnans in range(len(arrnans)):
                                                                              
        isnanarr = np.isnan(arrnans[qdnans])
        for qdnans2 in range(len(isnanarr)):   
            
            if(isnanarr[qdnans2]==True):
                
                newindex = np.arange(0,len(arrnans),1)
                
                for qdnans3 in range(len(newindex)):
                    arrnans[qdnans3][qdnans2] = -1e10
        
        newarrnanslist.append(arrnans[qdnans])

    for pindnan in range(len(newarrnanslist)):    
        newarrnanslist[pindnan] = np.array(newarrnanslist[pindnan])
        newarrnanslist[pindnan] = newarrnanslist[pindnan]\
                                  [newarrnanslist[pindnan]>-1e9]
    newarrnanslist = np.array(newarrnanslist)
    return newarrnanslist

#Ignore bad time intervals (BTIs)
def ignore_btis(arrays,tS,tE):
        
    tref = arrays[-1]
    oldarrays_list,newarrays_list = [[],[]]
                
    for qd in range(len(arrays)):
                                                                                
        #Identify BTIs
        for qd2 in range(len(tS)-1):
            
            btiS = tE[qd2]
            btiE = tS[qd2+1]
                        
            for qd3 in range(len(tref)):
                
                if(tref[qd3]>=btiS and tref[qd3]<=btiE\
                   and qd!=len(arrays)-1):
                    
                    newindarray = np.arange(0,len(arrays),1)
                    for qd4 in range(len(newindarray)):
                        arrays[qd4][qd3] = -1e10
                        
        #Remove BTIs
        if(qd!=len(arrays)-1):
            newarray = arrays[qd][arrays[qd]>-1e9]
            newarrays_list.append(newarray)
        oldarrays_list.append(arrays[qd])
        
    for pind in range(len(newarrays_list)):    
        newarrays_list[pind] = np.array(newarrays_list[pind])
    # newarrays_list = np.array(newarrays_list)

    for pind2 in range(len(oldarrays_list)):
        oldarrays_list[pind2] = np.array(oldarrays_list[pind2])
    # oldarrays_list = np.array(oldarrays_list)
        
    return oldarrays_list, newarrays_list
    
#Estimate covariance spectrum in time domain (Wilkinson and Uttley 2009)
def covariance_time_domain(lc,lcerr,reflc,reflcerr,Msegs):
        
    sigcov,sigxs_x,sigxs_y,sigerr_x,sigerr_y = [[],[],[],[],[]]
    Numpt = 0
    for arr in range(Msegs):
        
        #Split LC into M equal segments
        divs = int(len(reflc)/Msegs) 
        lctemp = lc[arr*divs:(arr+1)*divs]
        lctemperr = lcerr[arr*divs:(arr+1)*divs]
        reflctemp = reflc[arr*divs:(arr+1)*divs]
        reflctemperr = reflcerr[arr*divs:(arr+1)*divs]
        
        Numpt = len(lctemp)
        mutemp = np.mean(lctemp)
        mureftemp = np.mean(reflctemp)
                
        w1 = np.zeros(len(lctemp))
        w2 = np.zeros(len(lctemp))
        w3 = np.zeros(len(lctemp))
        w4 = np.zeros(len(lctemp))
        w5 = np.zeros(len(lctemp))

        for jp in range(len(lctemp)):
            w1[jp] = (lctemp[jp] - mutemp)*(reflctemp[jp] - mureftemp)
            w2[jp] = lctemperr[jp]**2
            w3[jp] = reflctemperr[jp]**2
            w4[jp] = (lctemp[jp]-mutemp)*(lctemp[jp]-mutemp)
            w5[jp] = (reflctemp[jp] - mureftemp)*(reflctemp[jp] - mureftemp)
                
        sigcov.append(np.mean(w1))
        sigerr_x.append(np.mean(w2))
        sigerr_y.append(np.mean(w3))
        sigxs_x.append(np.mean(w4) - np.mean(w2))
        sigxs_y.append(np.mean(w5) - np.mean(w3))
                
    mu_sigcov = np.mean(sigcov)
    mu_sigxs_x = np.mean(sigxs_x)
    mu_sigxs_y = np.mean(sigxs_y)
    mu_sigerr_x = np.mean(sigerr_x)
    mu_sigerr_y = np.mean(sigerr_y)
                
    mean_covariance = mu_sigcov/np.sqrt(mu_sigxs_y)
    cov_error = np.sqrt((mu_sigxs_x*mu_sigerr_y + mu_sigxs_y*mu_sigerr_x +\
    mu_sigerr_x*mu_sigerr_y)/(Msegs*Numpt*mu_sigxs_y))
                    
    if(np.isnan(mean_covariance)==True):
        mean_covariance = 0
        cov_error = 0
                
    return mean_covariance, cov_error

#Function to average PSD and CPSD over M segments
def Pbin(MsegPbin,freqsarr,Pxarr,Pyarr,Cxyarr):
        
    favg,Pxavg,dPxavg,Pyavg,dPyavg,Cxyavg,dCxyavg = [[],[],[],[],[],[],[]]
                            
    #Average CPSD and PSD over M segments
    for ipbin in range(np.shape(Cxyarr)[1]):
        
        px = 0
        py = 0
        cxy = 0
        freq = freqsarr[0][ipbin]
        
        for jpbin in range(np.shape(Cxyarr)[0]):
            
            cxy += Cxyarr[jpbin][ipbin]/MsegPbin
            px += Pxarr[jpbin][ipbin]/MsegPbin
            py += Pyarr[jpbin][ipbin]/MsegPbin
                
        favg.append(freq)
        Cxyavg.append(cxy)
        Pxavg.append(px)
        Pyavg.append(py)
        dPxavg.append(px/np.sqrt(MsegPbin))
        dPyavg.append(py/np.sqrt(MsegPbin))
        dCxyavg.append(cxy/np.sqrt(MsegPbin))
    
    favg = np.array(favg)
    Pxavg = np.array(np.real(Pxavg))
    Pyavg = np.array(np.real(Pyavg))
    dPxavg = np.array(np.real(dPxavg))
    dPyavg = np.array(np.real(dPyavg))
    Cxyavg = np.array(Cxyavg)
    dCxyavg = np.array(dCxyavg)
    
    #Return averaged quantities
    return favg,Pxavg,dPxavg,Pyavg,dPyavg,Cxyavg,dCxyavg

#Function to geometrically bin PSD and CPSD in frequency space
def fbin(bfact,farr,PXarr,PYarr,CXYarr,dPXarr,dPYarr,dCXYarr):
    
    avgf,avgPx,davgPx,avgPy,davgPy,avgCxy,davgCxy = [[],[],[],[],[],[],[]]
    Karr = []
        
    bmin = 0
    bmax = 0
    while(bmax<(len(farr))):
        
        if((bmax-bmin)!=0):
            bmin = bmax
        fmax = bfact*farr[bmax]        
        for index in range(bmin,len(farr)):
            if(farr[index]>=fmax):
                bmax = index
                break
            if(index==len(farr)-1 and farr[index]<=fmax):
                bmax = index + 1
                break
        if((bmax-bmin)==0):
            bfact+=0.1
                
        #Append relevant quantities        
        af = 0
        apx = 0
        apy = 0
        acxy = 0
        errpx = 0
        errpy = 0
        errcxy = 0
                
        for k5 in range(bmin,bmax):
            af += farr[k5]
            apx += PXarr[k5]
            apy += PYarr[k5]
            acxy += CXYarr[k5]
            errpx += dPXarr[k5]**2
            errpy += dPYarr[k5]**2
            errcxy += dCXYarr[k5]**2
        
        af/=(bmax-bmin)
        apx/=(bmax-bmin)
        apy/=(bmax-bmin)
        acxy/=(bmax-bmin)
        errpx = np.sqrt(errpx)/(bmax-bmin)
        errpy = np.sqrt(errpy)/(bmax-bmin)
        errcxy = np.sqrt(errcxy)/(bmax-bmin)
        
        avgf.append(af)
        avgPx.append(apx)
        avgPy.append(apy)
        avgCxy.append(acxy)
        davgPx.append(errpx)
        davgPy.append(errpy)
        davgCxy.append(errcxy)        
        Karr.append(bmax-bmin)
    
    #Return binned quantities
    return avgf,avgPx,avgPy,avgCxy,davgPx,davgPy,davgCxy,Karr

#Covariance spectrum (stingray)
def covariance_spectrum_stingray(evfile,\
    bwidth,fbmin,fbmax,refemin,refemax,egrid,normalisation):
                                 
    eventsst_ref = EventList.read(evfile,"hea",\
    additional_columns=["DET_ID"])
    frq_interval = [fbmin,fbmax]
    rf_band = [refemin,refemax]
    telapse_ref = eventsst_ref.time[-1]-eventsst_ref.time[0] 
    Msegstref = 0

    if (telapse_ref > 100*ks):
        Msegstref = 14

    if (telapse_ref > 50*ks and telapse_ref <= 100*ks):
        Msegstref = 10

    if (telapse_ref > 25*ks and telapse_ref <= 50*ks):
        Msegstref = 8

    if (telapse_ref > 15*ks and telapse_ref <= 25*ks):
        Msegstref = 3

    if (telapse_ref <= 15*ks):
        Msegstref = 1
     
    segsize_ref = telapse_ref/Msegstref    
    covspec = CovarianceSpectrum(eventsst_ref,\
              freq_interval=frq_interval,segment_size=segsize_ref,\
              bin_time=bwidth,energy_spec=egrid,\
              norm=normalisation,ref_band=rf_band)
    
    covspecE = covspec.energy
    covspecspt = covspec.spectrum
    covspecspterr = covspec.spectrum_error
            
    isnanarrcov = np.isnan(covspecspt)
    covspecE = covspecE[isnanarrcov==False]
    covspecspt = covspecspt[isnanarrcov==False]
    covspecspterr = covspecspterr[isnanarrcov==False]
    isnanarrcov = np.isnan(covspecspterr)
    covspecE = covspecE[isnanarrcov==False]
    covspecspt = covspecspt[isnanarrcov==False]
    covspecspterr = covspecspterr[isnanarrcov==False]

    return covspecE,covspecspt,covspecspterr

#Estimate time lag and covariance spectrum in Fourier domain (Uttley et al. 2014)
def time_lag_func(lc,lcerr,reflc,reflcerr,lcbkg,refbkg,ywindow,\
                  Mseg,bfactor,dt,stat):
    
    freqs,Px,Py,Cxy = [[],[],[],[]]
    
    # Compute noise level depending on whether the counting statistics
    # are Poissonian or not 
    fnyq = 0.5*(dt**-1)
    Pnoise = 0
    Prefnoise = 0
    Msegnew = 0

    # Average power spectrum and cross spectrum over M segments
    for k3t in range(Mseg):
        
        # Split LC into M equal segments
        div = int(len(reflc)/Mseg)
        lctemp = lc[k3t*div:(k3t+1)*div]
        lctemperr = lcerr[k3t*div:(k3t+1)*div]
        reflctemp = reflc[k3t*div:(k3t+1)*div]
        reflctemperr = reflcerr[k3t*div:(k3t+1)*div]
        lcbkgtemp = lcbkg[k3t*div:(k3t+1)*div]
        refbkgtemp = refbkg[k3t*div:(k3t+1)*div]
        ywindowtemp = ywindow[k3t*div:(k3t+1)*div]
        
        if(stat=="Poissonian"):
            
            Pnoise += (2*(np.mean(lctemp) +\
                       np.mean(lcbkgtemp))/(np.mean(lctemp))**2)
            Prefnoise += (2*(np.mean(reflctemp) +\
                          np.mean(refbkgtemp))/(np.mean(reflctemp))**2)
        
        if(stat!="Poissonian"):
            
            errsq = np.sum(lctemperr**2)/len(lctemperr)
            errrefsq = np.sum(reflctemperr**2)/len(reflctemperr)
            
            Pnoise += errsq/(fnyq*(np.mean(lctemp))**2)
            Prefnoise += errrefsq/(fnyq*(np.mean(reflctemp))**2) 
                
        if(np.sum(lctemp)>-100):
            
            # FFT of comparison-band LC
            Xn = fft.fft(lctemp) 
            fxn = fft.fftfreq(len(lctemp),d=dt)
            
            # FFT of reference-band LC
            Yn = fft.fft(reflctemp) 
            fyn = fft.fftfreq(len(reflctemp),d=dt)
            
            # FFT of window function
            Wn = fft.fft(ywindowtemp)
            fwn = fft.fftfreq(len(ywindowtemp),d=dt)
                                    
            Xn = Xn[fxn>0]
            Yn = Yn[fyn>0]
            Wn = Wn[fwn>0]
            fxn = fxn[fxn>0]
            fyn = fyn[fyn>0]
            fwn = fwn[fwn>0]
                                    
            # Compute PSD and CPSD with 
            # rms-squared normalisation for each segment
            normpsdx = (2.0*dt)/((len(lctemp))*(np.mean(lctemp))**2)
            normpsdy = (2.0*dt)/((len(reflctemp))*(np.mean(reflctemp))**2)
            normcross = (2.0*dt)/((len(lctemp))*(np.mean(lctemp))*\
                        (np.mean(reflctemp)))
                    
            Psdx = normpsdx*((np.conj(Xn))*(Xn))
            Psdy = normpsdy*((np.conj(Yn))*(Yn))
            Crossxy = normcross*((np.conj(Xn))*(Yn))
            Npsdx = len(Psdx)
                                                                    
            # Append CPSD and PSD for each segment to 
            # pass to functions for averaging and binning
            
            if(len(Crossxy)>0 and len(Psdx)>0 and len(Psdy)>0):
                                
                freqs.append(fxn)
                Px.append(Psdx)
                Py.append(Psdy)
                Cxy.append(Crossxy)
                Msegnew += 1
                        
            if(len(Crossxy)<=0 or len(Psdx)<=0 or len(Psdy)<=0):
                                                
                freqs.append(np.zeros(Npsdx))
                Px.append(np.zeros(Npsdx))
                Py.append(np.zeros(Npsdx))
                Cxy.append(np.zeros(Npsdx))
        
    #Average Pnoise and Prefnoise
    Pnoise /= Msegnew
    Prefnoise /= Msegnew
                    
    Mseg = Msegnew
    freqs = np.array(freqs)
    Px = np.array(Px)
    Py = np.array(Py)
    Cxy = np.array(Cxy)
                                
    # Average PSDs and CPSDs over M segements
    favg,Pxavg,dPxavg,Pyavg,dPyavg,Cxyavg,dCxyavg =\
    Pbin(Mseg,freqs,Px,Py,Cxy)
    
    fb = favg
    avgPx = Pxavg
    avgPy = Pyavg
    avgCxy = Cxyavg
    avgPxerr = dPxavg
    avgPyerr = dPyavg
    avgCxyerr = dCxyavg
    Karr = np.ones(len(Pxavg))
    
    # Implement frequency dependent binning of averaged PSDs and CPSDs
    if(bfactor>1):
                
        fb,avgPx,avgPy,avgCxy,avgPxerr,avgPyerr,avgCxyerr,Karr =\
        fbin(bfactor,favg,Pxavg,Pyavg,Cxyavg,dPxavg,dPyavg,dCxyavg)
    
    avgPx = np.real(avgPx)
    avgPy = np.real(avgPy)
    avgPx = np.array(avgPx)
    avgPxerr = np.array(avgPxerr)
    avgPy = np.array(avgPy)
    avgPyerr = np.array(avgPyerr)
    avgCxy = np.array(avgCxy)
    avgCxyerr = np.array(avgCxyerr)
    Karr = np.array(Karr)
    fb = np.array(fb)
    freqx = fb
    dfreqx = freqx[1]-freqx[0]
                                                        
    # Averaged number of samples
    nsamples = Mseg*Karr
        
    # Noise level of CPSD amplitude
    nbias = ((avgPx-Pnoise)*Prefnoise + (avgPy-Prefnoise)*Pnoise +\
            (Pnoise*Prefnoise))/nsamples
    
    # Complex-valued CPSD amplitude
    Cxyamp = (np.real(avgCxy))**2 + (np.imag(avgCxy))**2 - nbias
    
    avcxyreal = np.real(avgCxy)
    avcxyrealerr = np.real(avgCxyerr)
    avcxyimag = np.imag(avgCxy)
    avcxyimagerr = np.imag(avgCxyerr)
    
    dnbias = ((avgPxerr**2)*(Prefnoise**2) +\
              (avgPyerr**2)*(Pnoise**2))/(nsamples**2)
        
    dCxyamp = 4*((avcxyreal**2)*(avcxyrealerr**2) +\
                 (avcxyimag**2)*(avcxyimagerr**2)) + dnbias
        
    # Raw Coherence
    coherence = Cxyamp/((avgPx)*(avgPy))
        
    # Statistical uncertainty on raw coherence
    dcoherence = ((2.0/(nsamples))**(0.5))*(1 - coherence**2)/\
                 (abs(coherence))
    
    coherence = np.sqrt(coherence)
    dcoherence = 0.5*(dcoherence)/(coherence)
    intcoherence = Cxyamp/((avgPx-Pnoise)*(avgPy-Prefnoise)) #Intrinsic
        
    # Uncertainty in intrinsic coherence (from Vaughan and Nowak 1997)
    intcoherr = np.zeros(len(intcoherence))
    arbfact = 3 
            
    for u in range(len(coherence)):
                    
        # High powers, high measured coherence         
        cond1 = (arbfact*Pnoise)/(np.sqrt(nsamples[u]))
        cond2 = (arbfact*Prefnoise)/(np.sqrt(nsamples[u]))
        cond3 = (arbfact*nbias[u])/((avgPx[u])*(avgPy[u]))
        
        if((avgPx[u]-Pnoise)>cond1 and (avgPy[u]-Prefnoise)>cond2\
            and coherence[u]>cond3):
                            
            intcoherr[u] = ((nsamples[u])**-0.5)*\
                           (np.sqrt((2*nsamples[u]*nbias[u]**2)/\
                           (Cxyamp[u] - nbias[u])**2 +\
                           ((Pnoise)/(avgPx[u]-Pnoise))**2 +\
                           ((Prefnoise)/(avgPy[u]-Prefnoise))**2 +\
                           (nsamples[u]*dcoherence[u]**2)/\
                           (intcoherence[u]**2)))
                                
            intcoherr[u] *= intcoherence[u]
            intcoherr[u] = abs(intcoherr[u])
            if(np.isnan(intcoherr[u])=='True'):
                intcoherr[u] = 0
        
        # High powers, low measured coherence 
        else:    
            intcoherr[u] = np.sqrt(Prefnoise**2/(avgPx[u]-Prefnoise)**2/\
            nsamples[u] + Pnoise**2/(avgPy[u]-Pnoise)**2/\
            nsamples[u] + (dcoherence[u]/intcoherence[u])**2)
            intcoherr[u] *= intcoherence[u]
            intcoherr[u] = abs(intcoherr[u])
            if(np.isnan(intcoherr[u])=='True'):
                intcoherr[u] = 0
            
    # Compute phase lag as a function of frequency between the 
    # two energy bands        
    Cxyimag = np.imag(avgCxy)
    Cxyreal = np.real(avgCxy)
    phaselag = np.zeros(len(Cxyimag))
    dphaselag = np.zeros(len(Cxyimag))
                                            
    for hp3 in range(len(phaselag)):
                    
        #Clockwise
        phaselag[hp3] = np.arctan(Cxyimag[hp3]/Cxyreal[hp3])
        div = phaselag[hp3]/np.pi
        
        #Ensure phase lag is confined to between -pi to pi
        if(div>1):
                        
            divs = str(div)
            divnum = int(divs.split(".")[0])
                                            
            #Should be between 0 to 1
            if((divnum-1)%3==0):
                div -= divnum
            
            #Should be between 0 to 1
            if(divnum%3==0):
                div -= divnum

            #Should be between -1 to 0
            if((divnum+1)%3==0):
                div += (divnum+1)
            
        if(div<-1):
            
            divs = str(div)
            divnum = int(divs.split(".")[0])
            
            #Should be between 0 to 1
            if((divnum-1)%3==0):
                div -= (divnum-1)
            
            #Should be between 0 to 1
            if(divnum%3==0):
                div -= divnum

            #Should be between -1 to 0
            if((divnum+1)%3==0):
                div -= divnum
        
        coherence[hp3] = np.sqrt(coherence[hp3])
        phaselag[hp3] = div*np.pi
        dphaselag[hp3] = np.sqrt((1.0-coherence[hp3]**2)/\
                         (2.0*nsamples[hp3]*coherence[hp3]**2))
            
    return freqx,dfreqx,phaselag,dphaselag,coherence,dcoherence

#Estimate covariance spectrum in Fourier domain (Uttley et al. 2014)
def covariance_spectrum(tlc,lc,lcerr,reflc,reflcerr,lcbkg,refbkg,Mseg,\
                        bfactor,Dt,stats,fbmin,fbmax,window,rmbt):
                        
    #Compute noise level depending on whether the counting statistics
    #are Poissonian or not 
    fnyq = 0.5*(Dt**-1)
    Pnoise = 0
    Prefnoise = 0
    Msegnew = 0
    freqs,Pxn,Pyn,Cxyn,dPxn,dPyn,dCxyn = [[],[],[],[],[],[],[]]
        
    if(Mseg%2==0 and len(reflc)%2!=0):
                
        reflc = reflc[0:-1]
        reflcerr = reflcerr[0:-1]
        refbkg = refbkg[0:-1]
        
    #Average power spectrum and cross spectrum over M segments
    for kz in range(Mseg):
        
        #Split LC into M equal segments
        div = int(len(reflc)/Mseg)
        lctemp = lc[kz*div:(kz+1)*div]
        lctemperr = lcerr[kz*div:(kz+1)*div]
        reflctemp = reflc[kz*div:(kz+1)*div]
        reflctemperr = reflcerr[kz*div:(kz+1)*div]
        lcbkgtemp = lcbkg[kz*div:(kz+1)*div]
        refbkgtemp = refbkg[kz*div:(kz+1)*div]
        windowtemp = window[kz*div:(kz+1)*div]
                        
        #Ambient noise level in PSD
        if(stats=="Poissonian"):
            
            Pnoise += (2*(np.mean(lctemp) + np.mean(lcbkgtemp))/\
                      (np.mean(lctemp))**2)
            Prefnoise += (2*(np.mean(reflctemp) + np.mean(refbkgtemp))/\
                         (np.mean(reflctemp))**2)

        if(stats!="Poissonian"):
            
            errsq = 0
            errrefsq = 0
            
            for l in range(len(lctemperr)):
                errsq += lctemperr[l]**2
                errrefsq += reflctemperr[l]**2
                
            errsq/=len(lctemperr)
            errrefsq/=len(reflctemperr)
            
            Pnoise += errsq/(fnyq*(np.mean(lctemp))**2)
            Prefnoise += errrefsq/(fnyq*(np.mean(reflctemp))**2) 
        
        if(np.sum(lctemp)>0):
            
            lctemp = np.array(lctemp)
            lctemperr = np.array(lctemperr)
                                                            
            #FFT of comparison-band LC
            Xn = 0.5*(fft.fft(lctemp+lctemperr)+fft.fft(lctemp-lctemperr)) 
            Xnerr = 0.5*(fft.fft(lctemp+lctemperr)-fft.fft(lctemp-lctemperr))
            Xnconj = 0.5*(np.conj(Xn+Xnerr)+np.conj(Xn-Xnerr))
            Xnconjerr = 0.5*(np.conj(Xn+Xnerr)-np.conj(Xn-Xnerr))
            fxn = fft.fftfreq(len(lctemp),d=Dt)
                        
            #FFT of reference-band LC
            Yn = 0.5*(fft.fft(reflctemp+reflctemperr)+\
                 fft.fft(reflctemp-reflctemperr))
            Ynerr = 0.5*(fft.fft(reflctemp+reflctemperr)-\
                    fft.fft(reflctemp-reflctemperr))
            Ynconj = 0.5*(np.conj(Yn+Ynerr)+np.conj(Yn-Ynerr))
            Ynconjerr = 0.5*(np.conj(Yn+Ynerr)-np.conj(Yn-Ynerr))
            fyn = fft.fftfreq(len(reflctemp),d=Dt)
            
            # FFT of window function
            Wn = fft.fft(windowtemp)
            Wnconj = np.conj(Wn)
            
            if(rmbt=="True"):
                
                #Remove beats due to window
                Xn /= Wn
                Yn /= Wn
                Xnconj /= Wnconj
                Ynconj /= Wnconj
                        
            Xn = Xn[fxn>0] 
            Xnerr = Xnerr[fxn>0]
            Xnconj = Xnconj[fxn>0]
            Xnconjerr = Xnconjerr[fxn>0]
            Yn = Yn[fyn>0]
            Ynerr = Ynerr[fyn>0]
            Ynconj = Ynconj[fyn>0]
            Ynconjerr = Ynconjerr[fyn>0]
            Wn = Wn[fxn>0]
            Wnconj = Wnconj[fxn>0]
            fxn = fxn[fxn>0]
            fyn = fyn[fyn>0]
                                                                                                    
            # #Compute PSD and CPSD with 
            # rms-squared normalisation for each segment
            normpsdx = (2.0*Dt)/((len(lctemp))*(np.mean(lctemp))**2)
            normpsdw = (2.0*Dt)/((len(windowtemp))*(np.mean(windowtemp))**2)
            normpsdy = (2.0*Dt)/((len(reflctemp))*(np.mean(reflctemp))**2)
            normcross = (2.0*Dt)/((len(lctemp))*(np.mean(lctemp))*\
                        (np.mean(reflctemp)))
            
            #PSD
            Psdw = normpsdw*Wnconj*Wn
            Psdx = normpsdx*Xnconj*Xn
            dPsdx = normpsdx*(Xnconjerr*Xn + Xnconj*Xnerr)
            Psdy = normpsdy*Ynconj*Yn
            dPsdy = normpsdy*(Ynconjerr*Yn + Ynconj*Ynerr)
            Crossxy = normcross*Ynconj*Xn
            dCrossxy = normcross*(Ynconjerr*Xn + Ynconj*Xnerr)
                                    
            if(len(Crossxy)>0 and len(Psdx)>0 and len(Psdy)>0):
                                
                # Append CPSD and PSDs for each segment to 
                # pass to functions for averaging and binning
                freqs.append(fxn)
                Pxn.append(Psdx)
                Pyn.append(Psdy)
                Cxyn.append(Crossxy)
                Msegnew += 1
    
    #Average ref band and comp band noise powers
    Mseg = Msegnew
    Pnoise /= Mseg
    Prefnoise /= Mseg

    freqs = np.array(freqs)
    Pxn = np.array(Pxn)
    dPxn = np.array(dPxn)
    Pyn = np.array(Pyn)
    dPyn = np.array(dPyn)
    Cxyn = np.array(Cxyn)
    dCxyn = np.array(dCxyn)
                
    mean_covariance = 0
    err_mcov = 0

    if(len(Cxyn)>0):
                        
        # Average PSDs and CPSD over M segements
        freqx,Pxavg,dPxavg,Pyavg,dPyavg,Cxyavg,dCxyavg =\
        Pbin(Mseg,freqs,Pxn,Pyn,Cxyn)
        
        avgPx = Pxavg
        avgPy = Pyavg
        avgCxy = Cxyavg
        avgPxerr = dPxavg
        avgPyerr = dPyavg
        avgCxyerr = dCxyavg
        dfreqx = freqx[1]-freqx[0]
        Karr = np.ones(len(Pxavg))
        
        Karr = np.array(Karr)
        avgPxerr = np.array(avgPxerr)
        avgPyerr = np.array(avgPyerr)
        avgCxyerr = np.array(avgCxyerr)
        avgPx = np.array(avgPx)
        avgPy = np.array(avgPy)
        avgCxy = np.array(avgCxy)
        freqx = np.array(freqx)
        
        lgfreqx = np.log(freqx)
        lgavgPx = np.log(avgPx)
        lgavgPy = np.log(avgPy)
        lgavgPxerr = abs(avgPxerr/avgPx)
        lgavgPyerr = abs(avgPyerr/avgPy)

        # Averaged number of samples
        nsamples = Mseg*Karr
                                    
        # Noise level of CPSD
        nbias = ((avgPx-Pnoise)*Prefnoise + (avgPy-Prefnoise)*Pnoise +\
                (Pnoise*Prefnoise))/nsamples
                            
        # Compute CPSD amplitude from complex-valued cross spectrum
        Cxyamp = (np.real(avgCxy))**2 + (np.imag(avgCxy))**2 - nbias
        avcxyreal = np.real(avgCxy)
        avcxyrealerr = np.real(avgCxyerr)
        avcxyimag = np.imag(avgCxy)
        avcxyimagerr = np.imag(avgCxyerr)
        
        dnbias = ((avgPxerr**2)*(Prefnoise**2) +\
                  (avgPyerr**2)*(Pnoise**2))/(nsamples**2)
            
        Cxyamperr = 4*((avcxyreal**2)*(avcxyrealerr**2) +\
                    (avcxyimag**2)*(avcxyimagerr**2)) + dnbias
            
        # Raw Coherence
        coherence = Cxyamp/((avgPx)*(avgPy))
        intcoherence = Cxyamp/((avgPx-Pnoise)*(avgPy-Prefnoise))
        
        # Statistical uncertainty on raw coherence
        dcoherence = ((2.0/(nsamples))**(0.5))*(1 - intcoherence**2)/\
                     (abs(intcoherence))
        coherence = np.sqrt(coherence)
        dcoherence = 0.5*(dcoherence)/(coherence)
        
        # Compute covariance and its error over a desired frequency range
        rmsx = np.sqrt((avgPx-Pnoise)*(dfreqx)*(np.mean(lc))**2)
        rmsy = np.sqrt((avgPy-Prefnoise)*(dfreqx)*(np.mean(reflc))**2)
        rmsx_noise = np.sqrt((Pnoise)*((np.mean(lc))**2)*(dfreqx))
        rmsy_noise = np.sqrt((Prefnoise)*((np.mean(reflc))**2)*(dfreqx))
        covariance = (np.mean(lc))*np.sqrt((Cxyamp*dfreqx)/(avgPy-Prefnoise))
        dcovsqterm = ((covariance**2)*(rmsy_noise**2)+\
                      (rmsy**2)*(rmsx_noise**2)+\
                      (rmsy_noise**2)*(rmsx_noise**2))/(2*nsamples*rmsy**2)
        covariance_err = np.sqrt(dcovsqterm)

        # Filter over a specified frequency range
        covariance = covariance[freqx>fbmin]
        covariance_err = covariance_err[freqx>fbmin]
        coherence = coherence[freqx>fbmin]
        dcoherence = dcoherence[freqx>fbmin]
        avgCxy = avgCxy[freqx>fbmin]
        avgCxyerr = avgCxyerr[freqx>fbmin]
        freqxfilt = freqx[freqx>fbmin]
        covariance = covariance[freqxfilt<fbmax]
        covariance_err = covariance_err[freqxfilt<fbmax]
        coherence = coherence[freqxfilt<fbmax]
        dcoherence = dcoherence[freqxfilt<fbmax]
        avgCxy = avgCxy[freqxfilt<fbmax]
        avgCxyerr = avgCxyerr[freqxfilt<fbmax]
        freqxfilt = freqxfilt[freqxfilt<fbmax]
                        
        # Remove NANs
        isnancov = np.isnan(covariance)
        covariance = covariance[isnancov==False]
        coherence = coherence[isnancov==False]
        dcoherence = dcoherence[isnancov==False]
        covariance_err = covariance_err[isnancov==False]
        freqxfilt = freqxfilt[isnancov==False]
        
        isnancov = np.isnan(covariance_err)
        covariance = covariance[isnancov==False]
        covariance_err = covariance_err[isnancov==False]
        coherence = coherence[isnancov==False]
        dcoherence = dcoherence[isnancov==False]
        freqxfilt = freqxfilt[isnancov==False]
                                
        # Average over the above frequency range and propagate errors 
        # on covariance 
        mean_covariance = np.mean(covariance)
        err_mcov = 0
        for index in range(len(covariance)):
            err_mcov += covariance_err[index]**2
        err_mcov = np.sqrt(err_mcov)/len(covariance)
                
        if(len(covariance)==0):
            mean_covariance = 0
            err_mcov = 0
        
    return mean_covariance, err_mcov
               
#Cross spectrum (stingray)
def cross_spectrum_stingray(eventsst_ref,eventsst_comp,\
                            Msegstref,bwidth,normalisation,fbinmin,fbinmax):
    
    telapse_ref = eventsst_ref.time[-1]-eventsst_ref.time[0]
    segsize_ref = telapse_ref/Msegstref
    crossspec = AveragedCrossspectrum(eventsst_ref,eventsst_comp,\
                                      segment_size=segsize_ref,\
                                      norm=normalisation,dt=bwidth)
    cspec = crossspec.power
    cspecerr = crossspec.power_err
    cspecfreq = crossspec.freq
            
    cspec = cspec[cspecfreq>fbinmin]
    cspecerr = cspecerr[cspecfreq>fbinmin]
    cspecfreq = cspecfreq[cspecfreq>fbinmin]
    cspec = cspec[cspecfreq<fbinmax]
    cspecerr = cspecerr[cspecfreq<fbinmax]
    cspecfreq = cspecfreq[cspecfreq<fbinmax]
    
    cspecamp = cspec
    cspecamperr = cspecerr
                
    return cspecfreq,cspecamp,cspecamperr

#Generate LC from event-list
def make_lc(tarr,bintime,tstart,tstop,statlc):
    
    tobs = np.arange(tstart,tstop+bintime,bintime)
    robs = np.zeros(len(tobs))
    errobs = np.zeros(len(tobs))
    
    for k2 in range(len(tarr)):
        diffindex = np.argmin(abs(tarr[k2]-tobs))
        robs[diffindex] += 1
    
    if(statlc=='gauss'):
        errobs = np.std(robs)
    if(statlc=='poissonian'):
        errobs = np.sqrt(robs)
    
    robs/=bintime
    errobs/=bintime
    
    return tobs,robs,errobs

#Fill gaps in LC
def fgaps(arrsarrgap,constrategap,cthreshpoisson,binsize):
    
    reflc,ereflc,complc,ecomplc,wcomplc,tcomplc = arrsarrgap   
    
    ereflcpos = ereflc[reflc>0]
    complcpos = complc[reflc>0]
    ecomplcpos = ecomplc[reflc>0]
    tcomplcpos = tcomplc[reflc>0]
    reflcpos = reflc[reflc>0]
    Nsamples = len(reflcpos)
            
    timerefsim,reflccombsim,complccombsim,\
    errreflccombsim,errcomplccombsim = [[],[],[],[],[]]
    
    #Bootstrapped LC
    for jdsamp in range(Nsamples):
                                                        
        #Reference band and comparison band
        if(constrategap=="True"):
            
            mureflcconst = np.mean(reflcpos)
            mucomplcconst = np.mean(complcpos)
            
            muerrreflcconst = np.sqrt(np.sum(ereflcpos**2))/\
            len(reflcpos)
            muerrcomplcconst = np.sqrt(np.sum(ecomplcpos**2))/\
            len(complcpos)
                                    
            timerefsim.append(tcomplcpos[jdsamp])
            reflccombsim.append(mureflcconst)
            errreflccombsim.append(muerrreflcconst)
            complccombsim.append(mucomplcconst)
            errcomplccombsim.append(muerrcomplcconst)
        
        #Reference band
        randintref = np.random.randint(0,len(reflcpos),1)[0]
        ctsbinref = int(reflcpos[randintref]*binsize)

        if(constrategap=="False" and ctsbinref<cthreshpoisson):
           
            mureflc = np.random.poisson(ctsbinref,1)[0]/binsize 
            timerefsim.append(tcomplcpos[randintref])
            reflccombsim.append(mureflc)
            errreflccombsim.append(ereflcpos[randintref])

        if(constrategap=="False" and ctsbinref>=cthreshpoisson):
          
            mureflc = np.random.normal(reflcpos[randintref],\
                      ereflcpos[randintref],1)[0]
            timerefsim.append(tcomplcpos[randintref])
            reflccombsim.append(mureflc)
            errreflccombsim.append(ereflcpos[randintref])

        #Comparison band
        randintcomp = np.random.randint(0,len(complcpos),1)[0]
        ctsbincomp = int(complcpos[randintcomp]*binsize)
        
        if(constrategap=="False" and ctsbincomp<cthreshpoisson):
           
            mucomplc = np.random.poisson(ctsbincomp,1)[0]/binsize
            complccombsim.append(mucomplc)
            errcomplccombsim.append(ecomplcpos[randintcomp])

        if(constrategap=="False" and ctsbincomp>=cthreshpoisson):
           
            mucomplc = np.random.normal(complcpos[randintcomp],\
                       ecomplcpos[randintcomp],1)[0]
            complccombsim.append(mucomplc)
            errcomplccombsim.append(ecomplcpos[randintcomp])            
    
    timerefsim = np.array(timerefsim)
    reflccombsim = np.array(reflccombsim)
    errreflccombsim = np.array(errreflccombsim)
    complccombsim = np.array(complccombsim)
    errcomplccombsim = np.array(errcomplccombsim)    
    resultks = stats.ks_2samp(reflccombsim,reflc)
        
    timerefgap,reflccombgap,complccombgap,\
    errreflccombgap,errcomplccombgap,windowgap = [[],[],[],[],[],[]]
                        
    for jwinref in range(len(tcomplc)):
                                                                
        if(wcomplc[jwinref]==0):
            
            randnum = np.random.randint(0,len(tcomplcpos),1)[0]
                                    
            reflc[jwinref] = reflccombsim[randnum]
            ereflc[jwinref] = errreflccombsim[randnum]
            complc[jwinref] = complccombsim[randnum]
            ecomplc[jwinref] = errcomplccombsim[randnum]
            
            timerefgap.append(tcomplc[jwinref])
            reflccombgap.append(reflccombsim[randnum])
            errreflccombgap.append(errreflccombsim[randnum])
            complccombgap.append(complccombsim[randnum])
            errcomplccombgap.append(errcomplccombsim[randnum])
                              
    reflc = list(reflc)
    ereflc = list(ereflc)
    complc = list(complc)
    ecomplc = list(ecomplc)
    tcomplc,reflc,ereflc,\
    complc,ecomplc,wcomplc =\
    zip(*sorted(zip(tcomplc,reflc,ereflc,\
    complc,ecomplc,wcomplc)))
        
    timerefgap,reflccombgap,errreflccombgap,complccombgap,\
    errcomplccombgap = zip(*sorted(zip(timerefgap,reflccombgap,\
    errreflccombgap,complccombgap,errcomplccombgap)))
        
    tcomplc = np.array(tcomplc)
    reflc = np.array(reflc)
    ereflc = np.array(ereflc)
    complc = np.array(complc)
    ecomplc = np.array(ecomplc)
    wcomplc = np.array(wcomplc)
    timerefgap = np.array(timerefgap)
    reflccombgap = np.array(reflccombgap)
    errreflccombgap = np.array(errreflccombgap)
    complccombgap = np.array(complccombgap)
    errcomplccombgap = np.array(errcomplccombgap)
        
    return tcomplc, reflc, ereflc,\
    complc, ecomplc, wcomplc, timerefgap, reflccombgap, errreflccombgap,\
    complccombgap, errcomplccombgap

########################### Covariance spectrum ###############################

args = arguments()

# PSD parameters 
bfactor,reblog = args['geombin']
bfactor = float(bfactor)
reblog = float(reblog)
freqmin = args['freqmin']
freqmax = args['freqmax']
freqmin = float(freqmin)
freqmax = float(freqmax)
plotpsd = args['plotpsd']
statpow = args['statspsd']
splitscheme = args['splitscheme']
normpsd = args['normpsd']
gencov = args['covspec']
psdmods = args['powspecmod']

fminb = [freqmin]
fmaxb = [freqmax]
fminb = np.array(fminb)
fmaxb = np.array(fmaxb)
tthresh = (np.min(fminb))**-1

#LC parameters
plotlc = args['plotlc']
fillgaps,fillmethod = args['fillgaps']
segmentlc,tmin,tmax = args['segmentlc']
tmin *= ks
tmax *= ks

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

#MCMC simulations
runmcmc,Ntrialmcmc = args['runmcmc']

#Plot lags
plotlags = args['plotlags']

#Group covariance spectrum
groupscale = args['groupscale']
storagedir = "lags"

loc = os.getcwd()
os.chdir(loc + "/" + storagedir + "/")

keyobsid = "epn*net*obs*0*_1_*en3*comp*.lc"
obsidnum = []
for fobsid in sorted(glob.glob(keyobsid)):
    obsid = fobsid.split(".lc")[0].split("_")[2].split("obs")[1]
    obsidnum.append(obsid)
obsidnum = np.array(obsidnum)
obsidnum = np.unique(obsidnum)

Mseg = 1
ctsthreshpoisson = 10
siglag = 1.0
srcname = args['sourcename']
instarr = ["epn"]
labinst = ["EPN"]
col = ["bo","go","ro"]

comparecpsd = "False"
plotmcmc = "False"
plotcov = "False"
removebt = "False"
metmcmc = "timelags"

for kn in range(len(obsidnum)):   
                            
    Nenergies = 0
    for jn in sorted(glob.glob("epn*net*" + str(obsidnum[kn])+\
                               "*_1_*ref*.lc")):
        Nenergies += 1  
                    
    keyobs1 = "epn_net_obs*"
    keyobs2 = "*_1_*en3*ref.lc"
            
    for tempreflcfile in sorted(glob.glob(keyobs1+str(obsidnum[kn])+keyobs2)):
                                                                
        for ln in range(len(fminb)):
                                    
            cov,dcov,covtd,dcovtd,fvarc,fvarcerr = [[],[],[],[],[],[]]
    
            for qinstr in range(len(instarr)):
                
                instr = instarr[qinstr]
                enlag,denlag,mlag,mlagerr,mlagS,mlagerrS,\
                mufakelag,fakelagerr,energiesref,denergiesref,\
                = [[],[],[],[],[],[],[],[],[],[]]
            
                for k in range(Nenergies):
                                    
                    visnum = int(tempreflcfile.split("_")[3])
                    ennum = str(k+1)
                    ObsId = tempreflcfile.split("obs")[1].split("_")[0]
                    
                    refevfile = "epn_net_obs" + ObsId + "_" +\
                    str(visnum) + "_" + "en" + str(ennum) + "_ref.lc"
                    compevfile = "epn_net_obs" + str(ObsId) + "_" +\
                    str(visnum) + "_" + "en" + str(ennum) + "_comp.lc"
                    
                    #Reference band
                    reflcfile = "epn_net_obs" + str(ObsId) +\
                    "_" + str(visnum) + "_" + "en" + str(ennum) + "_ref.lc"
                                        
                    hdulistref = fits.open(reflcfile)
                    dataref = hdulistref[1].data
                    timeref = dataref['TIME']  
                    rateref = dataref['RATE']
                    errorref = dataref['ERROR']
                    tbeginR = hdulistref[2].header['TSTART']
                    tstartR = hdulistref[2].data['START']
                    tstopR = hdulistref[2].data['STOP']
                    bsizeref = timeref[1]-timeref[0]
                    telapse = timeref[-1]-timeref[0]
                                                                                
                    #Reference band (background)
                    refbkgfile = "epn_bkg_obs" + ObsId +\
                    "_" + str(visnum) + "_" + "en" + str(ennum) + "_ref.lc"
                    hdulistref_bkg = fits.open(refbkgfile)
                    dataref_bkg = hdulistref_bkg[1].data
                    timerefbkg = dataref_bkg['TIME']
                    raterefbkg = dataref_bkg['RATE']
                    errorrefbkg = dataref_bkg['ERROR']
                    
                    #Comparison band
                    complcfile = "epn_net_obs" + str(ObsId) + "_" +\
                    str(visnum) + "_" + "en" + str(ennum) + "_comp.lc"
                    hducomp = fits.open(complcfile)  
                    timecomp = hducomp[1].data['TIME']
                    ratecomp = hducomp[1].data['RATE']
                    errorcomp = hducomp[1].data['ERROR']
                    obsidcomp = hducomp[0].header['OBS_ID']
                                                        
                    #Comparison band (background)         
                    compbkgfile = "epn_net_obs" + str(ObsId) + "_" +\
                    str(visnum) + "_" + "en" + str(ennum) + "_comp.lc" 
                    hdubkgcomp = fits.open(compbkgfile)
                    ratecompbkg = hdubkgcomp[1].data['RATE']
                    errorcompbkg = hdubkgcomp[1].data['ERROR']
                    tstartC = hducomp[2].data['START']
                    tstopC = hducomp[2].data['STOP']
                    
                    #Additional information
                    raobj = hdulistref_bkg[0].header['RA_OBJ']
                    decobj = hdulistref_bkg[0].header['DEC_OBJ']
                    telescope = hdulistref_bkg[0].header['TELESCOP']
                    filterobs = hdulistref_bkg[0].header['FILTER']
                    inst = hdulistref_bkg[0].header['INSTRUME']
                    dmobs = hdulistref_bkg[0].header['DATE-OBS'].split("T")[0]
                    dmend = hdulistref_bkg[0].header['DATE-END'].split("T")[0]
                    tmobs = hdulistref_bkg[0].header['DATE-OBS'].split("T")[1]
                    tmend = hdulistref_bkg[0].header['DATE-END'].split("T")[1]
                    
                    #Remove NANs
                    arraysR =\
                    np.transpose(np.column_stack((rateref,\
                    errorref,ratecomp,errorcomp,raterefbkg,\
                    errorrefbkg,ratecompbkg,errorcompbkg,\
                    timeref,timecomp)))
                    arraysR = remove_nans_lc(arraysR)
                    rateref,\
                    errorref,ratecomp,errorcomp,raterefbkg,\
                    errorrefbkg,ratecompbkg,errorcompbkg,\
                    timeref,timecomp = arraysR
                
                    #Add rectangular window (reference-band)                  
                    arraysW = np.transpose(np.column_stack((timeref,rateref)))
                    windowref = rect_window(arraysW,tstartR,tstopR)
                    
                    #Add rectangular window (comparison-band)                  
                    arraysWc = np.transpose(np.column_stack((timeref,\
                                                             ratecomp)))
                    windowcomp = rect_window(arraysWc,tstartC,tstopC)
                                                            
                    infilecov = "covflux" + str(ln+1) +\
                    "_" + str(ObsId) + ".dat"
                    outfilecov = "covspec" + str(ln+1) +\
                    "_" + str(ObsId) + ".pha"
                    groupfilecov = "covspec_grouped_" + str(ln+1) + "_" +\
                    str(ObsId) + ".pha"
                                                                                                                                                                                    
                    #Choose Mseg depending on exposure time
                    if(splitscheme=="True"):
                                                                                    
                        if (telapse > 100*ks):
                            Mseg = 15
                    
                        if (telapse > 50*ks and telapse <= 100*ks):
                            Mseg = 12
                    
                        if (telapse > 25*ks and telapse <= 50*ks):
                            Mseg = 8
                    
                        if (telapse > 15*ks and telapse <= 25*ks):
                            Mseg = 3
                    
                        if (telapse <= 15*ks):
                            Mseg = 1
                                        
                    quantcomp = hducomp[1].header['DSVAL6']
                    
                    if(quantcomp!='TABLE'):
                                                                        
                        energymin = float(quantcomp.split(",")[0]\
                                       .split(":")[0])/1000.0
                        energymax = float(quantcomp.split(",")[0]\
                                       .split(":")[1])/1000.0
                        energiesmean = 0.5*(energymax+energymin)
                        denergiesmean = 0.5*(energymax-energymin)
                        
                        energiesref.append(energiesmean)
                        denergiesref.append(denergiesmean)
                        enlag.append(energiesmean)
                        denlag.append(denergiesmean)            
                    
                    if(quantcomp=='TABLE'):
                                                
                        quantcomp = hducomp[1].header['DSVAL5']
                                                                        
                        energymin = float(quantcomp.split(",")[0]\
                                       .split(":")[0])/1000.0
                        energymax = float(quantcomp.split(",")[0]\
                                       .split(":")[1])/1000.0
                        energiesmean = 0.5*(energymin+energymax)
                        denergiesmean = 0.5*(energymin-energymax)
                        
                        energiesref.append(energiesmean)
                        denergiesref.append(denergiesmean)
                        enlag.append(energiesmean)
                        denlag.append(denergiesmean)
                                                                                                                                                                
                    if(len(rateref)>0):
                        
                        windowcomb = np.array(windowref)
                        reflccomb = np.array(rateref)
                        errreflccomb = np.array(errorref)
                        reflcbkgcomb = np.array(raterefbkg)
                        errreflcbkgcomb = np.array(errorrefbkg)
                        complccomb = np.array(ratecomp)
                        errcomplccomb = np.array(errorcomp)
                        complcbkgcomb = np.array(ratecompbkg)
                        errcomplcbkgcomb = np.array(errorcompbkg)
                        timecombref = bsizeref*np.arange(0,len(reflccomb),1)
                                                                        
                        #Multiply by window function
                        reflccomb *= windowcomb
                        complccomb *= windowcomb

                        #Compare CPSDs
                        if(comparecpsd=="True"):
                            
                            #Stingray LCs and events
                            countcomp = complccomb*bsizeref
                            cerrorcomp = errcomplccomb*bsizeref
                            countref = reflccomb*bsizeref
                            cerrorref = errreflccomb*bsizeref
                            countcompbkg = complcbkgcomb*bsizeref
                            countrefbkg = raterefbkg*bsizeref
                            LcComp = Lightcurve(timecomp,countcomp,\
                            error=cerrorcomp,dt=bsizeref)
                            LcRef = Lightcurve(timecomp,countref,\
                            error=cerrorref,dt=bsizeref)
                            evcomp = EventList.from_lc(LcComp)   
                            evref = EventList.from_lc(LcRef)
                                                        
                            #CPSD (stingray)
                            cspecfrq, crosspec, crosspecerr =\
                            cross_spectrum_stingray(evref,evcomp,\
                            Mseg,bsizeref,normpsd,fminb,fmaxb)
                    
                            #CPSD
                            fxy, cpxy, cpxyerr =\
                            covariance_spectrum(complccomb,\
                            errcomplccomb,reflccomb,errreflccomb,\
                            complcbkgcomb,reflcbkgcomb,\
                            Mseg,bfactor,bsizeref,statpow,fminb,fmaxb,\
                            windowref)
                                                
                            plt.figure()
                            plt.subplot(211)
                            plt.errorbar(fxy,np.real(cpxy),\
                                         yerr=abs(np.real(cpxyerr)),\
                                         fmt='g.')
                            plt.errorbar(cspecfrq,np.real(crosspec),\
                                         yerr=abs(np.real(crosspecerr)),\
                                         fmt='m.')
                            plt.subplot(212)
                            plt.errorbar(fxy,np.imag(cpxy),\
                                         yerr=abs(np.imag(cpxyerr)),\
                                         fmt='g.')
                            plt.errorbar(cspecfrq,np.imag(crosspec),\
                                         yerr=abs(np.imag(crosspecerr)),\
                                         fmt='m.')
                            plt.show()
                
                        if(segmentlc=="True"):
                            
                            complccomb = complccomb[timecombref>=tmin]
                            errcomplccomb = errcomplccomb[timecombref>=tmin]
                            reflccomb = reflccomb[timecombref>=tmin]
                            errreflccomb = errreflccomb[timecombref>=tmin]
                            complcbkgcomb = complcbkgcomb[timecombref>=tmin]
                            errcomplcbkgcomb =\
                            errcomplcbkgcomb[timecombref>=tmin]
                            reflcbkgcomb = reflcbkgcomb[timecombref>=tmin]
                            errreflcbkgcomb =\
                            errreflcbkgcomb[timecombref>=tmin]
                            windowcomb = windowcomb[timecombref>=tmin]
                            timecombref = timecombref[timecombref>=tmin]
                            
                            complccomb = complccomb[timecombref<=tmax]
                            errcomplccomb = errcomplccomb[timecombref<=tmax]
                            reflccomb = reflccomb[timecombref<=tmax]
                            errreflccomb = errreflccomb[timecombref<=tmax]
                            complcbkgcomb = complcbkgcomb[timecombref<=tmax]
                            errcomplcbkgcomb =\
                            errcomplcbkgcomb[timecombref<=tmax]
                            reflcbkgcomb = reflcbkgcomb[timecombref<=tmax]
                            errreflcbkgcomb =\
                            errreflcbkgcomb[timecombref<=tmax]
                            windowcomb = windowcomb[timecombref<=tmax]
                            
                            timecombref =\
                            bsizeref*(np.arange(0,len(windowcomb),1))
                                                                        
                        if(psdmods=="True"):
                                                                      
                            #Generate a GPR-based PSD model                                
                            lgfreqypsd,lgavgPypsd,lgavgPyerrpsd,\
                            lgfreqypsdmod,lgymodpsd,mulcpsd =\
                            psdmodgp(timecombref,reflccomb,errreflccomb,\
                            complccomb,errcomplccomb,reflcbkgcomb,\
                            complcbkgcomb,Mseg,bfactor,bsizeref,\
                            statpow,removebt)
    
                            # plt.errorbar(lgfreqypsd,lgavgPypsd,\
                            #              yerr=lgavgPyerrpsd,fmt='k.')
                            # plt.plot(lgfreqypsdmod,lgymodpsd,'b-')
                            # plt.show()
                                                                                                                                                                                                                                           
                        if(fillgaps=="True"):
                            
                            if(timmerkoenig=='True'):
                                
                                #Generate a GPR-based PSD model
                                #(Reference-band)
                                lgfreqypsd,lgavgPypsd,lgavgPyerrpsd,\
                                lgfreqypsdmod,lgmodpsdy,mulcpsdy =\
                                psdmodgp(timecombref,reflccomb,errreflccomb,\
                                complccomb,errcomplccomb,reflcbkgcomb,\
                                complcbkgcomb,Mseg,bfactor,bsizeref,\
                                statpow,removebt)
                                    
                                #Generate a GPR-based PSD model
                                #(Comparison-band)
                                lgfreqxpsd,lgavgPxpsd,lgavgPxerrpsd,\
                                lgfreqxpsdmod,lgmodpsdx,mulcpsdx =\
                                psdmodgp(timecombref,complccomb,\
                                errcomplccomb,reflccomb,\
                                errreflccomb,complcbkgcomb,\
                                reflcbkgcomb,Mseg,bfactor,bsizeref,\
                                statpow,removebt)
                                                                            
                                if(len(lgfreqypsd)>0):
                                    
                                    freqypsd = np.exp(lgfreqypsd)
                                    psdmody = 10**lgmodpsdy
                                    psdmodx = 10**lgmodpsdx
                                
                                    telapsecomb =\
                                    timecombref[-1]-timecombref[0]
                                    
                                    freqgenmin = np.min(freqypsd)
                                    freqgenmax = np.max(freqypsd)
                                                                                
                                    timesimref,refsimlc,errrefsimlc =\
                                    genlc(bsizeref,telapsecomb,freqgenmin,\
                                    freqgenmax,psdmody,mulcpsdy)
                                                                            
                                    timesimref,compsimlc,errcompsimlc =\
                                    genlc(bsizeref,telapsecomb,freqgenmin,\
                                    freqgenmax,psdmodx,mulcpsdx)
                                                                            
                            if(bootstrap=="True"):
                                                                                      
                                #Add a floor
                                if(np.min(reflccomb)<0):
                                    reflccomb += abs(np.min(reflccomb))
                                if(np.min(complccomb)<0):
                                    complccomb += abs(np.min(complccomb))
                                timecombref -= timecombref[0]
                                
                                #Remove NANs
                                arraysR =\
                                np.transpose(np.column_stack((reflccomb,\
                                errreflccomb,complccomb,\
                                errcomplccomb,windowcomb,timecombref)))
                                
                                timecombref,reflccomb,errreflccomb,\
                                complccomb,errcomplccomb,windowcomb,\
                                timesimref,refsimlc,errrefsimlc,compsimlc,\
                                errcompsimlc =\
                                fgaps(arraysR,stdwin,ctsthreshpoisson,\
                                      bsizeref)

                        if(fillgaps=="False"):
                            
                            #Remove BTIs and NANs from reference band and 
                            #comparison band
                            arraysR =\
                            np.transpose(np.column_stack((reflccomb,\
                            errreflccomb,reflcbkgcomb,errreflcbkgcomb,\
                            complccomb,errcomplccomb,\
                            complcbkgcomb,errcomplcbkgcomb,timecombref)))
                            arraysR = remove_nans_lc(arraysR)
                                                            
                            #Ignore gaps
                            if(len(tstartR)>1):
                                                                                
                                arraysR, arraysN =\
                                ignore_btis(arraysR,tstartR,tstopR)
                                reflccomb,errreflccomb,\
                                reflcbkgcomb,errreflcbkgcomb,\
                                complccomb,errcomplccomb,\
                                complcbkgcomb,errcomplcbkgcomb = arraysN
                            
                            if(len(tstartR)==1):
                                                                
                                reflccomb,errreflccomb,\
                                reflcbkgcomb,errreflcbkgcomb,\
                                complccomb,errcomplccomb,\
                                complcbkgcomb,errcomplcbkgcomb,\
                                timecombref = arraysR
    
                        if(plotlc=="True" and ln==0):
                                                                                    
                            #Plot LCs
                            laben = str(float(energymin)) + "-" +\
                            str(float(energymax))
                            labelsrccomp = "Comparison band LC: " +\
                            laben + " keV"    
                                                            
                            plt.figure()
                            plt.title("XMM-Newton (EPIC-PN) lightcurves")
                            plt.errorbar(timecombref/ks,reflccomb,\
                                         yerr=errreflccomb,\
                                         fmt='r.')
                            plt.errorbar(timecombref/ks,complccomb,\
                                         yerr=errcomplccomb,\
                                         label=labelsrccomp,\
                                         fmt='b.')
                            
                            
                            if(fillgaps=="True" and fillmethod=="B"):
                                
                                plt.errorbar(timesimref/ks,refsimlc,\
                                yerr=errrefsimlc,fmt='g.',\
                                label="Interpolated: Bootstrapped")
                                plt.errorbar(timesimref/ks,compsimlc,\
                                yerr=errcompsimlc,fmt='k.',\
                                label="Interpolated: Bootstrapped")
                                    
                            plt.plot(timecombref/ks,windowcomb,'m-')
                            plt.tick_params(axis='both', which='major',\
                                            labelsize=14)
                            plt.legend(loc="upper right")
                            plt.ylabel("Count rate [s$^{-1}$]",\
                                       fontsize=14)
                            plt.xlabel("Time [ks]",fontsize=14)
                            plt.show()
    
                            #Mean and RMS of rate
                            muref = np.mean(rateref)
                            dmuref = np.sum(errorref**2)/len(rateref)
                            mucomp = np.mean(ratecomp)
                            dmucomp = np.sum(errorcomp**2)/len(ratecomp)
                            
                        mucomp = np.mean(complccomb)
                        muref = np.mean(reflccomb)
                        countcomp = complccomb*bsizeref
                        countcomp_err = errcomplccomb*bsizeref
                        countref = reflccomb*bsizeref
                        countref_err = errreflccomb*bsizeref
                        countcompbkg = ratecompbkg*bsizeref
                        countrefbkg = raterefbkg*bsizeref
                        
                        #Add a floor
                        countcomp -= np.min(countcomp)
                        countref -= np.min(countref)
                        
                        timecomp = np.arange(0,len(countcomp),1)*bsizeref
                        timeref = timecomp
                        ywindowR = np.ones(len(rateref))
                                                                        
                        #Compute event lists from LC
                        lccomp = Lightcurve(timecombref,countcomp,\
                        error=countcomp_err,dt=bsizeref)
                        lcref = Lightcurve(timeref,countref,\
                        error=countref_err,dt=bsizeref)
                        
                        segsizest =\
                        (lcref.time[-1]-lcref.time[0]+bsizeref)/Mseg
                        
                        evcomplc = EventList.from_lc(lccomp)
                        evreflc = EventList.from_lc(lcref)
                                                
                        #Compute time lags using Stingray
                        csa = AveragedCrossspectrum.from_events(evcomplc,\
                        evreflc, segment_size=segsizest,\
                        dt=bsizeref, norm="abs",use_common_mean=True,\
                        silent=True)
                            
                        rebloglags = 0.0
                        csa = csa.rebin_log(rebloglags)
                        csaamp = csa.power
                        csaamperr = csa.power_err
                        freq_lag = csa.freq
                        coh, coh_e = csa.coherence()
                        lag, lag_e = csa.time_lag()
                        lagp, lag_ep = csa.phase_lag()
    
                        csumref = bsizeref*np.sum(lcref.countrate)
                        csumreferr = bsizeref*\
                        np.sqrt(np.sum(lcref.countrate_err**2))
                        
                        psdcomp =\
                        AveragedPowerspectrum.from_lightcurve(lccomp,\
                        segment_size=segsizest,norm="frac",silent=True)
                        
                        psdref =\
                        AveragedPowerspectrum.from_lightcurve(lcref,\
                        segment_size=segsizest,norm="frac",silent=True)
                        
                        psdref = psdref.rebin_log(reblog)
                        psdcomp = psdcomp.rebin_log(reblog)
                
                        # Model PSDs
                        
                        # Initialise fitting parameters
                        pfitref = psdref.power
                        perrfitref = psdref.power_err
                        freqfitref = psdref.freq
                        pfit = psdcomp.power
                        perrfit = psdcomp.power_err
                        freqfit = psdcomp.freq
                    
                        # Perform PSD fitting in log space and transform back
                        pfitreflog = np.log10(pfitref)
                        perrfitreflog = (perrfitref)/(pfitref*np.log(10.0))
                        ampinitref = np.mean(pfitreflog)
                        alphainitref = 2.0
                        ampbesinitref = 100.0
                        norderinitref = 1
                        paramsinitial =\
                        [ampinitref,alphainitref,ampbesinitref,norderinitref]
                        fitobj = kmpfit.Fitter(residuals=resid_besselmod,\
                        data=(freqfitref,pfitreflog,perrfitreflog))
                        fitobj.fit(params0=paramsinitial)
                        chi2min = fitobj.chi2_min
                        dof = fitobj.dof
                        
                        ampbestref = fitobj.params[0]
                        ampbestreferr = fitobj.xerror[0]
                        alphabestref = fitobj.params[1]
                        alphabestreferr = fitobj.xerror[1]
                        ampbesbestref = fitobj.params[2]
                        ampbesbestreferr = fitobj.xerror[2]
                        norderbestref = fitobj.params[3]
                        norderbestreferr = fitobj.xerror[3]

                        bestfitparsref =\
                        [ampbestref,alphabestref,ampbesbestref,norderbestref]
                                       
                        fmod = np.linspace(np.min(psdref.freq),\
                        np.max(psdref.freq),1000)
                        pmodref = 10**(besselmod(bestfitparsref,fmod))
                                            
                        pfitlog = np.log10(pfit)
                        perrfitlog = (perrfit)/(pfit*np.log(10.0))
                        ampinitcomp = np.mean(pfitlog)
                        alphainitcomp = 2.0
                        ampbesinitcomp = 1.0
                        norderinitcomp = 1
                        paramsinitialcomp = [ampinitcomp,alphainitcomp,\
                        ampbesinitcomp,norderinitcomp]
                        fitobjcomp = kmpfit.Fitter(residuals=resid_besselmod,\
                        data=(freqfit,pfitlog,perrfitlog))
                        fitobjcomp.fit(params0=paramsinitialcomp)
                        chi2mincomp = fitobjcomp.chi2_min
                        dofcomp = fitobjcomp.dof
                        
                        ampbestcomp = fitobjcomp.params[0]
                        ampbestcomperr = fitobjcomp.xerror[0]
                        alphabestcomp = fitobjcomp.params[1]
                        alphabestcomperr = fitobjcomp.xerror[1]
                        ampbesbestcomp = fitobjcomp.params[2]
                        ampbesbestcomperr = fitobjcomp.xerror[2]
                        norderbestcomp = fitobjcomp.params[3]
                        norderbestcomperr = fitobj.xerror[3]

                        bestfitparscomp = [ampbestcomp,alphabestcomp,\
                                           ampbesbestcomp,norderbestcomp]
                        pmodref = 10**(besselmod(bestfitparsref,fmod))
                        pmodcomp = 10**(besselmod(bestfitparscomp,fmod))
                                                
                        #Plot power-spectra
                        if(plotpsd=="True"):
                            
                            plt.figure()
                            plt.title("PSD " + str(srcname) +\
                                      " Obs ID: " + str(ObsId))
                            plt.errorbar(psdref.freq,psdref.power,\
                                         yerr=psdref.power_err,fmt='r.',\
                                         label="Reference band")
                            # plt.plot(fmod,pmodref,'k--',\
                            #          label="PL fit [reference band]")
                            plt.errorbar(psdcomp.freq,psdcomp.power,\
                                         yerr=psdcomp.power_err,fmt='b.',\
                                         label="Comparison band")
                            # plt.plot(fmod,pmodcomp,'k-',\
                            #          label="PL fit [comparison band]")
                            plt.yscale("log")
                            plt.xscale("log")
                            plt.xlabel("Frequency [Hz]",fontsize=14)
                            plt.ylabel("Power (fractional rms) [Hz$^{-1}$]",\
                                       fontsize=14)
                            plt.tick_params(axis='both', which='major',\
                                            labelsize=14)
                            plt.legend(loc="best")
                            plt.show()
    
                        #Fractional variability
                        Fracvarref = 0.5*(integrate.simpson(psdref.power+\
                                          psdref.power_err,\
                                          psdref.freq) +\
                                          integrate.simpson(psdref.power-\
                                          psdref.power_err,\
                                          psdref.freq))
                        dFracvarref = 0.5*(integrate.simpson(psdref.power+\
                                           psdref.power_err,\
                                           psdref.freq) -\
                                           integrate.simpson(psdref.power-\
                                           psdref.power_err,\
                                           psdref.freq))
                
                        Fracvarcomp = 0.5*(integrate.simpson(psdcomp.power+\
                                           psdcomp.power_err,\
                                           psdcomp.freq) +\
                                           integrate.simpson(psdcomp.power-\
                                           psdcomp.power_err,\
                                           psdcomp.freq))
                        dFracvarcomp = 0.5*(integrate.simpson(psdcomp.power+\
                                            psdcomp.power_err,\
                                            psdcomp.freq) -\
                                            integrate.simpson(psdcomp.power-\
                                            psdcomp.power_err,\
                                            psdcomp.freq))
                        Fracvarref = np.sqrt(Fracvarref)
                        dFracvarref = 0.5*(dFracvarref/Fracvarref)
                        Fracvarcomp = np.sqrt(Fracvarcomp)
                        dFracvarref = 0.5*(dFracvarcomp/Fracvarcomp)
    
                        # Compute time lags (lag-energy spectrum) 
                        # using my method
                        bfactorlags = 1.0
                        freqS, dfreqS, lagS, lag_eS, cohS, coh_eS =\
                        time_lag_func(reflccomb,errreflccomb,\
                        complccomb,errcomplccomb,reflcbkgcomb,\
                        complcbkgcomb,windowcomb,\
                        Mseg,bfactorlags,bsizeref,stats)
                        
                                                                                                
                        # Compute time lags (lag-frequency spectrum)
                        # using my method
                        bfactorlags = bfactor
                        FreqS, dFreqS, lagfreqS, lagfreq_eS,\
                        cohfreqS, cohfreq_eS =\
                        time_lag_func(reflccomb,errreflccomb,\
                        complccomb,errcomplccomb,reflcbkgcomb,\
                        complcbkgcomb,windowcomb,\
                        Mseg,bfactorlags,bsizeref,stats)
                                                
                        lagfreqS = lagfreqS/(2.0*np.pi*FreqS)
                        lagfreq_eS = lagfreq_eS/(2.0*np.pi*FreqS)
                        
                        if(plotlags=="True"):
                        
                            plt.errorbar(FreqS,-lagfreqS,yerr=lagfreq_eS,\
                                         fmt='ko')
                            plt.errorbar(FreqS,-lagfreqS,yerr=lagfreq_eS,\
                            markersize=8,marker='o',linestyle='dotted')                            
                            plt.xscale("log")
                            plt.xlim(7e-5,2e-3)
                            plt.xlabel("Frequency [Hz]",\
                                       fontsize=14)
                            plt.title("Lag frequency spectrum [" +\
                                      str(srcname) + "]",fontsize=14)
                            plt.ylabel("Time lag [s]",fontsize=14)
                            plt.tick_params(axis='both',\
                            which='major',labelsize=16)
                            plt.show()
                    
                        # Compute covariance
                        # Time domain
                        intcovtd,intcoverrtd =\
                        covariance_time_domain(complccomb,errcomplccomb,\
                                               reflccomb,errreflccomb,\
                                               Mseg)
                        covtd.append(intcovtd)
                        dcovtd.append(intcoverrtd)
                        
                        # Frequency domain                
                        intcov,intcoverr =\
                        covariance_spectrum(\
                        timecombref,complccomb,errcomplccomb,\
                        reflccomb,errreflccomb,complcbkgcomb,\
                        reflcbkgcomb,Mseg,bfactor,bsizeref,\
                        statpow,fminb[ln],fmaxb[ln],windowcomb,removebt)
                
                        cov.append(intcov)
                        dcov.append(intcoverr)
                        
                        if(runmcmc=="True"):
                            
                            fakelags = []
                            
                            for z0 in range(int(Ntrialmcmc)):
                                
                                tlagfk =\
                                mcmc_det(bsizeref,telapse,Mseg,bfactor,\
                                fminb[ln],fmaxb[ln],\
                                ampbestref,alphabestref,ampbesbestref,\
                                norderbestref,ampbestcomp,alphabestcomp,\
                                ampbesbestcomp,norderbestcomp,\
                                mucomp,muref,plotmcmc,statpow,metmcmc)
                                    
                                tlagfk = np.array(tlagfk)
                                fakelags.append(tlagfk)
                                                    
                            fakelags = np.array(fakelags)
                            mufakelag.append(np.mean(fakelags))
                            fakelagerr.append(np.std(fakelags))
                                                
                        if(len(lagS)>0):
                            
                            #Filter over a frequency range
                            lag = lag[freq_lag>=fminb[ln]]
                            lag_e = lag_e[freq_lag>=fminb[ln]]
                            lagp = lagp[freq_lag>=fminb[ln]]
                            lag_ep = lag_ep[freq_lag>=fminb[ln]]
                            freq_lag = freq_lag[freq_lag>=fminb[ln]]
                            lag = lag[freq_lag<=fmaxb[ln]]
                            lag_e = lag_e[freq_lag<=fmaxb[ln]]
                            lagp = lagp[freq_lag<=fmaxb[ln]]
                            lag_ep = lag_ep[freq_lag<=fmaxb[ln]]
                            freq_lag = freq_lag[freq_lag<=fmaxb[ln]]
                            
                            lagS = lagS[freqS>=fminb[ln]]
                            lag_eS = lag_eS[freqS>=fminb[ln]]
                            cohS = cohS[freqS>=fminb[ln]]
                            coh_eS = coh_eS[freqS>=fminb[ln]]
                            freqS = freqS[freqS>=fminb[ln]]
                            lagS = lagS[freqS<=fmaxb[ln]]
                            lag_eS = lag_eS[freqS<=fmaxb[ln]]
                            cohS = cohS[freqS<=fmaxb[ln]]
                            coh_eS = coh_eS[freqS<=fmaxb[ln]]
                            freqS = freqS[freqS<=fmaxb[ln]]
                                                        
                            #Remove NANs                                        
                            arrayslag =\
                            np.transpose(np.column_stack((lagS,lag_eS,\
                            cohS,coh_eS,freqS,lag,lagp,lag_e,\
                            lag_ep,freq_lag)))
                            arrayslag = remove_nans_lags(arrayslag)
                            lagS,lag_eS,cohS,coh_eS,freqS,lag,\
                            lagp,lag_e,lag_ep,\
                            freq_lag = arrayslag
                                                                                                                                                                                                        
                        residlag =\
                        (lagS - lag)/(np.sqrt(lag_eS**2)+np.sqrt(lag_e**2))
                        residerr = np.ones(len(residlag))
                                                                    
                        if(len(lagS)==0):
                                                
                            mean_fbS = np.mean(freqS) 
                            mean_lagS = 0
                            mean_lagSerr = 0
                            mean_fb = np.mean(freq_lag)
                            mean_lag = 0
                            mean_lagerr = 0
                            
                            mlagS.append(np.nan)
                            mlagerrS.append(np.nan)
                            mlag.append(np.nan)
                            mlagerr.append(np.nan)
    
                        if(len(lagS)>0):
                                                                                                        
                            #Phase wrapping (stingray)
                            for pq in range(len(lagS)):
                                
                                #Shift by pi 
                                sigthresh = 1.2
                                ediff =\
                                np.sqrt(lag_ep[pq]**2 + lag_eS[pq]**2)
                                diff = (lagp[pq] - lagS[pq])/ediff 
                                niter = 5
                                kiter = 0
                                                
                                while(abs(diff)>sigthresh):
                                    
                                    diff = (lagp[pq] - lagS[pq])/ediff 
                                    
                                    if(diff<0 and abs(diff)>sigthresh):
                                        lagS[pq] -= np.pi
                                        
                                    if(diff>0 and abs(diff)>sigthresh):
                                        lagS[pq] += np.pi
                                    
                                    kiter += 1
                                    
                                    if(kiter>niter or abs(diff)<sigthresh):
                                        break
                                                            
                            tlagS = lagS/(2.0*np.pi*freqS)
                            tlag_eS = lag_eS/(2.0*np.pi*freqS)
                            tlag = lagp/(2.0*np.pi*freq_lag)
                            tlag_e = lag_ep/(2.0*np.pi*freq_lag)
                                                                                    
                            mean_fbS = np.mean(freqS) 
                            mean_lagS = np.median(tlagS)
                            mean_lagSerr =\
                            np.sqrt(np.sum(tlag_eS**2))/len(tlagS)
                            mean_fb = np.mean(freq_lag)
                            mean_lag = np.median(tlag)
                            mean_lagerr =\
                            np.sqrt(np.sum(tlag_e**2))/len(tlag)
                                                        
                            mlagS.append(mean_lagS)
                            mlagerrS.append(mean_lagSerr)
                            mlag.append(mean_lag)
                            mlagerr.append(mean_lagerr)
                                                                            
            enlag = np.array(enlag)
            mlagS = np.array(mlagS)
            mlag = np.array(mlag)
            denlag = np.array(denlag)
            mlagerrS = np.array(mlagerrS)
            mlagerr = np.array(mlagerr)  
            cov = np.array(cov)
            dcov = np.array(dcov)
            energiesref = np.array(energiesref)
            denergiesref = np.array(denergiesref)
                                                                        
            mufakelag = np.array(mufakelag)
            fakelagerr = np.array(fakelagerr)
            confint1 = mufakelag - siglag*(0.5*fakelagerr)
            confint2 = mufakelag + siglag*(0.5*fakelagerr)
                              
            if(plotlags=="True"):
                fnamesave = "lag_energy_" + str(ObsId) + ".dat"
                Z = np.column_stack((enlag,denlag,mlag/ks,mlagerr/ks))
                np.savetxt(fnamesave,Z,fmt='%s',delimiter='  ')
    
            lab1 = "Frequency band: (" + str(fminb[ln]) + "-" +\
                   str(fmaxb[ln]) + ") Hz, Instrument: " +\
                   str("EPIC-pn")
            ylim1 = (np.min(mlagS)-np.std(mlagS))/ks
            ylim2 = (np.max(mlagS)+np.std(mlagS))/ks
            
            if(np.sum(cov)>0 and np.sum(dcov)>0\
               and gencov=="True"):
                
                if(plotcov=="True"):
                                    
                    plt.figure()
                    plt.errorbar(energiesref,cov,yerr=dcov,fmt='b.',\
                                 label="Method 1")
                    plt.errorbar(energiesref,covtd,yerr=dcovtd,fmt='g.',\
                                 label="Method 2")
                    plt.legend(loc="best")
                    plt.show()
                
                #Unfold spectrum using instrumental response      
                rmffile = "epn_" + str(ObsId) +\
                          "_" + str(visnum) + ".rmf"
                ancrfile = "epn_" + str(ObsId) +\
                           "_" + str(visnum) + ".arf"
                specfile = "epn_spec" + str(visnum) +\
                           "_grp_" + str(ObsId) + ".fits"
        
                hdulist2 = fits.open(specfile)
                header2 = hdulist2[1].header
                backscal = header2['BACKSCAL']
                corrscal = header2['CORRSCAL']
                areascal = header2['AREASCAL']
                backfile = "NONE"
                
                # Convert to XSPEC readable format using response file
                hdulist = fits.open(rmffile)
                data = hdulist[2].data
                channel = data['CHANNEL']
                EMIN = data['E_MIN']
                EMAX = data['E_MAX']
                ndetchans = len(channel)
                chans = np.arange(1,ndetchans+1,1)
                energiesref = np.array(energiesref)
                denergiesref = np.array(denergiesref)
                fluxes = np.zeros(len(chans))
                dfluxes = np.zeros(len(chans))
                
                for j3 in range(len(cov)):
                    for k3p in range(len(chans)):
                        if(0.5*(EMIN[k3p]+EMIN[k3p])>=energiesref[j3]):
                            fluxes[chans[k3p]] = cov[j3]
                            dfluxes[chans[k3p]] = dcov[j3]
                            break
            
                Qcov = np.column_stack((chans,fluxes,dfluxes))
                np.savetxt(infilecov,Qcov,fmt='%i %s %s',delimiter='   ')
                                                
                comm_unfold = "ascii2pha infile=" + infilecov +\
                " outfile=" + outfilecov +\
                " chanpres=yes dtype=2 rows=- qerror=yes tlmin=1 detchans=" +\
                str(ndetchans) + " telescope=" + str(telescope) +\
                " instrume=" + str(inst) + " detnam=EPIC-PN" +\
                " filter=" + str(filterobs) + " phaversn=1.1.0 " +\
                "exposure=" + str(telapse/Mseg) + " backscal=" +\
                str(backscal) + " backfile=" + backfile +\
                " corrscal=" + str(corrscal) +\
                " corrfile=NONE areascal=" + str(areascal) +\
                " ancrfile=" + ancrfile + " respfile=" + rmffile +\
                " date_obs=" + str(dmobs) + " time_obs=" + str(tmobs) +\
                " date_end=" + str(dmend) + " time_end=" + str(tmend) +\
                " ra_obj=" + str(raobj) + " dec_obj=" + str(decobj) +\
                " equinox=2000.0 hduclas2=TOTAL chantype=PI clobber=yes"
                os.system(comm_unfold)
                    
                #Group spectrum using ftgrouppha
                comm_group = "ftgrouppha infile=" +\
                outfilecov + " backfile=" + backfile +\
                " respfile=" + rmffile +\
                " outfile=" + groupfilecov +\
                " grouptype=optsnmin groupscale=" +\
                str(groupscale) + " minchannel=-1 maxchannel=-1"
                os.system(comm_group)
                        
            if(np.sum(cov)>0 and np.sum(dcov)>0):
    
                #Lag-energy spectrum
                if(plotlags=="True"):
                                        
                    if(ln==0):
                        fig = plt.figure(figsize=(8,6))
                        
                    subplt = int(str(len(fminb)) + '1' + str(ln+1))
                    
                    ax1 = fig.add_subplot(subplt)
                    enlag = np.array(enlag)
                    mlagS = np.array(mlagS)
                    mlag = np.array(mlag)
                    denlag = np.array(denlag)
                    mlagerrS = np.array(mlagerrS)
                    mlagerr = np.array(mlagerr)
                    mufakelag = np.array(mufakelag)
                    fakelagerr = np.array(fakelagerr)
                                                                                                    
                    confint1 = mufakelag - 3*0.5*fakelagerr
                    confint2 = mufakelag + 3*0.5*fakelagerr
                    
                    lab1 = "Frequency band: (" + str(fminb[ln]) + "-" +\
                    str(fmaxb[ln]) + ") Hz, Instrument: " +\
                    str(labinst[qinstr])
                         
                    fname_save = "lag_energy_obsid_" + str(ObsId) + ".png"
                    ylim1 = (np.min(mlagS)-np.std(mlagS))/ks
                    ylim2 = (np.max(mlagS)+np.std(mlagS))/ks
                    fnamesave = "lag_energy" + str(ln+1) + ".dat"
                    
                    Z = np.column_stack((enlag,denlag,mlag/ks,mlagerr/ks))
                    np.savetxt(fnamesave,Z,fmt='%s',delimiter='  ')
                                                            
                    #Lag-energy spectrum
                    if(ln==0):
                        ax1.set_title("Lag-energy spectrum " +\
                        str(srcname) + ":" +\
                        " ObsID " + str(ObsId),fontsize=18)
                                        
                    ax1.errorbar(enlag,mlagS/ks,xerr=abs(denlag),\
                                 yerr=mlagerrS/ks,\
                    fmt=col[qinstr],alpha=0.5,label=lab1,\
                    markersize=8,marker='o',linestyle='dotted')
                    
                    # ax1.errorbar(enlag,mlag/ks,\
                    #              xerr=abs(denlag),yerr=mlagerr/ks,\
                    # fmt='k.',alpha=1.0,label="Stingray",\
                    # markersize=8,marker='o',linestyle='dotted')
                        
                    ax1.tick_params(axis='both', which='major',labelsize=18)
                    ax1.set_xscale("log")
                    ax1.get_xaxis().\
                    set_major_formatter(matplotlib.ticker.ScalarFormatter())
                    ax1.get_xaxis().get_major_formatter().labelOnlyBase =\
                    False  
                    ax1.set_xticks([0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
    
                    if(runmcmc=="True"):
                                                                                                                
                        ax1.plot(enlag,mufakelag/ks,'k--')
                        ax1.fill_between(enlag,confint1/ks,confint2/ks,\
                        alpha=0.25,\
                        label=r"1$\sigma$ confidence level [from MCMC]")
                        ax1.legend(loc="best",shadow=False,framealpha=0.2)
                                   
                    if(qinstr==0):
                        ax1.set_ylabel("Time lag [ks]",fontsize=18)
                        
                    ax1.set_xlabel("Energy [keV]",fontsize=18)
                    plt.savefig("lag_espec_" + str(ObsId) + ".png", dpi=100)
                    plt.subplots_adjust(hspace=0)
        plt.show()
        

    
