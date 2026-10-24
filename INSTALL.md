To clone the repository:
git clone https://github.com/raj1294/XcovSpec

Then navigate to the XcovSpec directory: N.B: The directory "must" be named XcovSpec (don't change this name):
cd XcovSpec

Run install.sh:

chmod u+x install.sh
./install.sh 

Note: Users must have pre-installed anaconda, XMM-SAS (https://www.cosmos.esa.int/web/xmm-newton/sas-installation), HEASOFT and XSPEC (with version compatible with the user’s system (see https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html). Try conda install —name xcovspec heasoft, however I would try to install from source via the official webpages. The developer does not take responsibility for installing the above packages. Please note that this version of XcovSpec is designed to be compatible with Mac Apple Darwin 24.4.0 M2, although earlier or later versions will work provided python v3.10 can be installed. Please don’t remove any contents from the install.sh file.

When you run install.sh, the terminal will install conda packages and request the user to proceed with [y/n] from time to time. Once the installation finishes successfully, you can run the code. Set the following PYTHONPATH before doing so

PYTHONPATH=$HOME/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}

Then simply edit the run.sh file:

python XcovSpec.py -srcname "NGC 5204 X-1" -obsids 0405690501 -refemin 0.3 -refemax 12.0 -mincts 30 -minbcts 30 -dtbinqpo 0.07336496 -dtbincov 100 -dtbinbkg 100 -srad 0.0083 -brad 0.0083 -srcdet 1e-16 -rmflares False -texp 5 -plc True -plags True -ppsd True -split False -statpower False -normpower abs -flgaps True,B -seglc True,0,500 -fmin 1e-4 -fmax 5e-4 -gbin 1.0,0.0 -gscale 1 -rmcmc False,1000 -gencov True -psdmods False -egrid minmax,0.3,10.0,4 -aflag True -psearch False -bsub False

Then finally:

chmod u+x run.sh
./run.sh


To uninstall the package: 

conda remove -n xcovspec --all, followed by conda clean -tp

