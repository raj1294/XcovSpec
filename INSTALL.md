To clone the repository:

git clone https://github.com/raj1294/XcovSpec

Then navigate to the XcovSpec directory: N.B: The directory "must" be named XcovSpec (don't change this name). Note that users must have pre-installed anaconda, XMM-SAS (https://www.cosmos.esa.int/web/xmm-newton/sas-installation), HEASOFT and XSPEC (with version compatible with the user’s system (see https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html). 

cd XcovSpec

chmod u+x install.sh

./install.sh 

Then, activate the xcovspec conda environment, and be sure to set the PYTHONPATH.

source activate xcovspec

PYTHONPATH=$HOME/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}

Edit the run.sh file and change parameters to suit your observation goals.

chmod u+x run.sh

python XcovSpec.py -srcname "NGC 1313 X-1" -obsids 0803990301 -refemin 0.3 -refemax 12.0 -mincts 50 -minbcts 50 -dtbinqpo 100 -dtbincov 100 -dtbinbkg 100 -srad 0.0083 -brad 0.0083 -srcdet 1e-16 -rmflares False -texp 5 -split False -statpower False -normpower abs -flgaps True,B -seglc True,0,500 -fmin 1e-4 -fmax 5e-4 -gbin 1.1,0.1 -gscale 1 -rmcmc False,1000 -plc False -plags False -ppsd False -gencov True -egrid minmax,0.3,10.0,15 -aflag True -psearch True -bsub True -stdir lag_freq

Run the code:

./run.sh

During the run, the terminal will open ds9 and ask the user to place the background region file(s) in the "proc" directory of the associated observation. If there are two visits, please save two region files, named bkg1.reg and bkg2.reg. If that goes correctly, the code should produce lag-energy spectra and covariance spectra, as directed by the user above. Remember that you can set the energy grid, reference band, choose any source and observation identifier from the XMM DR4 catalogue and run MCMC simulations to estimate the lag-significance using the Timmer and Koenig approach, extract light-curves for pulsation search etc. Finally, to uninstall the package: 

conda remove -n xcovspec --all, followed by conda clean -tp

N.B: The user must have pre-installed heasoft and headas (with xspec enabled) inside the xcovspec environment prior to running the code. Try:

conda install heasarc heasoft (inside the xcovspec environment)

I would try to install heasoft from source via the official webpages. The developer does not take responsibility for installing heasoft. Please note that this version of XcovSpec is designed to be compatible with Mac Apple Darwin 24.4.0 M2, although earlier or later versions will work provided python v3.10 can be installed. Please don’t remove contents from the install.sh file. When you run install.sh, the terminal will install conda packages and request the user to proceed with [y/n]. Once the installation finishes successfully, you can run the code. 

Also, be aware that one needs to initialise heasoft prior to running covariance.py, ensuring that HEADAS points to the correct directory. Also, "unset PFILES" as this ensures there is no conflict between HEADAS and CIAO. One can also run individual scripts from intermediate stages of the pipeline by commenting out the xcovspec.py file.
