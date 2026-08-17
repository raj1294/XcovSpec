conda create --name xcovspec python==3.10
conda install --name xcovspec -c https://cxc.cfa.harvard.edu/conda/ciao ciao caldb_main
#conda install —name XcovSpec heasoft
source activate xcovspec
python -m pip install -q astropy scipy scikit-learn stingray numpy==2.0.0 numba==0.67.0 cython wget astromy-ds9 astroquery
cd /Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages
wget https://www.astro.rug.nl/software/kapteyn/kapteyn-3.4.tar.gz
tar -xvf kapteyn-3.4.tar.gz
rm -f kapteyn-3.4.tar.gz
cd kapteyn-3.4
export CFLAGS="-Wno-error=int-conversion"
python setup.py install
cd ~/software/XcovSpec/
source deactivate
