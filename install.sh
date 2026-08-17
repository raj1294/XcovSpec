conda create --name xcovspec python==3.10
conda install --name xcovspec -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao caldb_main ciao-contrib sherpa ds9
source activate xcovspec
python -m pip install astropy scipy scikit-learn stingray numpy==2.0.0 numba==0.67.0 cython wget astroquery
cd /Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages
wget https://www.astro.rug.nl/software/kapteyn/kapteyn-3.4.tar.gz
tar -xvf kapteyn-3.4.tar.gz
rm -f kapteyn-3.4.tar.gz
cd kapteyn-3.4
export CFLAGS="-Wno-error=int-conversion"
python setup.py install
cd ~/software/XcovSpec/
python -m pip install numpy==2.0.0
PYTHONPATH=/Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}
