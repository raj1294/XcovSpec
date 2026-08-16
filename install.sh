conda create --name xcovspec python==3.10
conda install --name xcovspec -c https://cxc.cfa.harvard.edu/conda/ciao ciao caldb_main
#conda install —name XcovSpec heasoft
source activate xcovspec
python -m pip install astropy scipy scikit-learn stingray numpy==2.0.0 numba==0.67.0 cython wget astromy-ds9 astroquery
python -m pip install https://www.astro.rug.nl/software/kapteyn/kapteyn-3.4.tar.gz
export PYTHONPATH=/opt/homebrew/lib/python3.10/site-packages:${PYTHONPATH}
