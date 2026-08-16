Copyright 2026 Rajath Sathyaprakash 

XcovSpec version 1.0 (relased August 2026) 

Please add citation to GitHub repository after using this code.

The uploaded set of python-based command line programs are intended to quickly and reliably decipher meaningful data trends from the XMM-Newton archive via the extraction of light-curves, power density spectra, images, time-delay profiles and covariance spectra. These five data measures are often used to characterise variability in many different classes of sources including accreting X-ray binaries, AGN, novae and galaxy clusters.  

With the extraction of light-curves the user can control the time bin size and the energy grid either via the ‘array’ approach, in which case the energy grid is written explicitly separated by commas, or via the ‘minmax’ approach, in which case the user only need specify the minimum and maximum energy bin and the number of energy bins in logarithmic space. The scripts also enable the extraction of power-density spectra (PDS) with geometric binning enabled and tunable normalisation options (absolute-rms normalisation, fractional rms-squared normalisation or leahy normalisation), options to discard background flares from light-curves using the approach described in the XMM-Newton SAS tutorial threads, and the computation of lag-energy spectra, lag-frequency spectra, time-averaged spectra and covariance spectra, with the user able to pre-filter the light-curve over a narrow frequency band, interpolate across gaps due to bad-time intervals via the Timmer-Koenig approach or via bootstrapped light-curves. When computing lag-energy profiles between light-curves of two different energy bands for a real dataset, the user can choose to run MCMC simulations (simply a boolean switch) to assess the significance of the lags. Indeed, the user can control whether or not to remove flagged pixels when extracting filtered event files or time-averaged spectra. 

The tool automatically downloads the required data files from the XMM-Newton archive, based on the input source name. The source name needs to be a valid registry of the XMM-Newton observatory. When downloading the data, one can specify only the required Obs IDs in order to prevent the pipeline from downloading all the available datasets, which could potentially exceed the local hard-drive capacity. The tool runs a source detection algorithm in order to apply a spatial filter (mainly used to centroid the location of point-spread-function) to extract the required datasets by approximating the PSF with ELLBETA profiles. 

To install :

Simply run ./install.sh on your bash terminal.

Users must have pre-installed XMM-SAS from the following website: https://www.cosmos.esa.int/web/xmm-newton/sas-installation and the latest heasoft version including xspec (see https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html: maybe try conda install —name xcovspec heasoft, but try to install from source) compatible with the user’s system. Please note that this version of XcovSpec is designed to be compatible with {\textbf{Apple Darwin 24.4.0}}, although earlier versions will work provided python v3.10 can be installed. Please don’t remove any contents from the install.sh file, especially the PYTHONPATH directory.

Simply run the following python file as follows:

python ~/Documents/SEAWIND/code/Communication/XcovSpec.py -srcname “Ark 564”  -obsids 0861600101,0861600201 -refemin 0.3 -refemax 12.0 -mincts 100 -minbcts 100 -dtbinqpo 100 -dtbincov 100 -dtbinbkg 100 -srad 0.0083 -brad 0.0083 -srcdet 1e-16 -rmflares False -texp 5 -plc False -plags False -ppsd True -split False -statpower False -normpower frac -flgaps True,B -seglc True,0,500 -fmin 1e-4 -fmax 5e-4 -gbin 1.0,0.0 -gscale 1 -rmcmc False,1000 -gencov True -psdmods False -egrid minmax,0.3,10.0,16 -aflag True -psearch True -bsub True

In the above case, the source name is “Ark 564”, with the specific ObsIDs downloaded (i.e. 0861600101,0861600201), the reference energy band light-curve spanning the 0.3-12.0 keV energy range, with at-least 100 counts grouped for the time averaged spectrum (after subtracting the background), with the time bin size selected to be 100 s (for the pulsations and/or QPO search) and 100 s for the lag-energy spectra and covariance spectra, source and background extraction radius set to 30 arc-seconds, with the background flares not chosen to be excluded and a threshold of at-least 5 ks required to continue with the analysis. Please consult the help file for more detail on the specific command-line arguments.

usage: XcovSpec.py [-h] -srcname SOURCENAME -plc PLOTLC -plags PLOTLAGS -ppsd PLOTPSD -split SPLITSCHEME -statpower STATSPSD -normpower NORMPSD -refemin
                   REFERENCE_ENERGY_MIN -refemax REFERENCE_ENERGY_MAX -flgaps FILLGAPS -seglc SEGMENTLC -fmin FREQMIN -fmax FREQMAX -gbin GEOMBIN
                   -gscale GROUPSCALE -rmcmc RUNMCMC -psdmods POWSPECMOD -egrid ENERGY_GRID -gencov COVSPEC -mincts MINIMUM_CTS -minbcts MINIMUM_CTS_BKG
                   -srad SRCRAD -brad BKGRAD -rmflares REMOVE_BKG_FLARES -srcdet SIGTHRESH -dtbinbkg BKG_BIN_TIME -dtbincov BINNING_TIME_COV -aflag
                   ADD_FLAG -obsids OBSERVATION_IDS [-texp THRESHOLD_EXP_TIME] -dtbinqpo BINNING_TIME_QPO -psearch PULSE_SEARCH -bsub BKG_SUB

Generate covariance spectra, lag-energy spectra plot power-spectral densities and plot light-curves

options:

  -h, --help            show this help message and exit
  
  -srcname SOURCENAME, --sourcename SOURCENAME Target Name
  
  -plc PLOTLC, --plotlc PLOTLC Plot LC? [Enter either True or False]
  
  -plags PLOTLAGS, --plotlags PLOTLAGS Plot Lag-energy spectrum? [Enter either True or False]
  
  -ppsd PLOTPSD, --plotpsd PLOTPSD Plot PSD? [Enter either True or False]
  
  -split SPLITSCHEME, --splitscheme SPLITSCHEME Segment LC? [Enter either True or False]
  
  -statpower STATSPSD, --statspsd STATSPSD Does the LC follow a Poisson distribution? [Enter either True or False]
  
  -normpower NORMPSD, --normpsd NORMPSD Enter normalisation of PSD [abs or frac or leahy]
  
  -refemin REFERENCE_ENERGY_MIN, --reference_energy_min REFERENCE_ENERGY_MIN Reference-band minimum energy [keV]
  
  -refemax REFERENCE_ENERGY_MAX, --reference_energy_max REFERENCE_ENERGY_MAX Reference-band maximum energy [keV]
                        
  -flgaps FILLGAPS, --fillgaps FILLGAPS Fill LC gaps? [Enter 2 values separated by comma of type: Boolean(Enter either True or False) String(Interpolation scheme: Enter B (bootstrapping) , T (timmer-koenig) or S (window deconvolution))
  
  -seglc SEGMENTLC, --segmentlc SEGMENTLC Segment LC? [Enter 3 values separated by comma of type: True/False (Boolean), Start of LC in ks (float), End of LC in ks (float)]
                        
  -fmin FREQMIN, --freqmin FREQMIN Minimum Fourier Frequency for Covariance [Enter a floating point value in Hz]
                        
  -fmax FREQMAX, --freqmax FREQMAX Maximum Fourier Frequency for Covariance [Enter a floating point value in Hz]
                        
  -gbin GEOMBIN, --geombin GEOMBIN Geometric binning factor [Enter two floating point values]
                        
  -gscale GROUPSCALE, --groupscale GROUPSCALE Scaling factor to group Covariance spectrum [Enter an integer value]
                        
  -rmcmc RUNMCMC, --runmcmc RUNMCMC Run MCMC? [Enter 2 values separated by comma of type: Boolean(Enter either True or False), Int(Enter number of MCMC simulations)
  -psdmods POWSPECMOD, --powspecmod POWSPECMOD Model power spectral density? [Enter either True or False]
                        
  -egrid ENERGY_GRID, --energy_grid ENERGY_GRID Specify comparison-band energy grid [in keV]
                        
  -gencov COVSPEC, --covspec COVSPEC Generate Covariance Spectrum? [Enter either True or False]
                        
  -mincts MINIMUM_CTS, --minimum_cts MINIMUM_CTS Minimum counts in source spectrum
                        
  -minbcts MINIMUM_CTS_BKG, --minimum_cts_bkg MINIMUM_CTS_BKG Minimum counts in background spectrum
                        
  -srad SRCRAD, --srcrad SRCRAD Target Name
                        
  -brad BKGRAD, --bkgrad BKGRAD Background extraction radius [in degrees]
                        
  -rmflares REMOVE_BKG_FLARES, --remove_bkg_flares REMOVE_BKG_FLARES Remove background flares?
                        
  -srcdet SIGTHRESH, --sigthresh SIGTHRESH Source detection threshold
                        
  -dtbinbkg BKG_BIN_TIME, --bkg_bin_time BKG_BIN_TIME Background binning time [in seconds]
                        
  -dtbincov BINNING_TIME_COV, --binning_time_cov BINNING_TIME_COV Binning time of LC for covariance spectra [in seconds]
                        
  -aflag ADD_FLAG, --add_flag ADD_FLAG Add Flag to extract events?
                        
  -obsids OBSERVATION_IDS, --observation_ids OBSERVATION_IDS List of Observation Identifiers
                        
  -texp THRESHOLD_EXP_TIME, --threshold_exp_time THRESHOLD_EXP_TIME Only use observation if elapsed time > tthresh [in ks]
                        
  -dtbinqpo BINNING_TIME_QPO, --binning_time_qpo BINNING_TIME_QPO Binning time of LC for QPO searches [in seconds]
                        
  -psearch PULSE_SEARCH, --pulse_search PULSE_SEARCH Extract light-curves with DT=DTMIN=73ms?
                        
  -bsub BKG_SUB, --bkg_sub BKG_SUB Use epiclccorr to subtract background?

  Upcoming: Tools to model lag-energy spectra, lag-frequency spectra, time-averaged spectra and light-curves etc.
