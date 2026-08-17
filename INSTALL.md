To download: 

Clone repository from github using:
git clone https://github.com/raj1294/XcovSpec

To install :

cd to the XcovSpec directory and simply run the following commands on your bash terminal: 
chmod u+x install.sh
./install.sh 

N.B: Users must have pre-installed anaconda, XMM-SAS (https://www.cosmos.esa.int/web/xmm-newton/sas-installation) and heasoft (version compatible with the user’s system, including xspec: see https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/download.html). Try conda install —name xcovspec heasoft, however I would try to install from source via the official webpages. The developer does not take responsibility for installing the above packages. 

Please note that this version of XcovSpec is designed to be compatible with Apple Darwin 24.4.0, although earlier or later versions will work provided python v3.10 can be installed. Please don’t remove any contents from the install.sh file.

Set the following PYTHONPATH before running the code
PYTHONPATH=/Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}

To run the code: Simply run the following python file as follows from your terminal

python XcovSpec.py -srcname "NGC 5204 X-1" -obsids 0405690501 -refemin 0.3 -refemax 12.0 -mincts 30 -minbcts 30 -dtbinqpo 0.07336496 -dtbincov 100 -dtbinbkg 100 -srad 0.0083 -brad 0.0083 -srcdet 1e-16 -rmflares False -texp 5 -plc True -plags True -ppsd True -split False -statpower False -normpower abs -flgaps True,B -seglc True,0,500 -fmin 1e-4 -fmax 5e-4 -gbin 1.0,0.0 -gscale 1 -rmcmc False,1000 -gencov True -psdmods False -egrid minmax,0.3,10.0,4 -aflag True -psearch False -bsub False
