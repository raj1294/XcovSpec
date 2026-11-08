conda create --name xcovspec python==3.10
source activate xcovspec
conda install --name xcovspec -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao caldb_main ciao-contrib sherpa mpmath llvmlite numba
python -m pip install astropy scipy scikit-learn stingray==1.1.2.4 numpy==1.26.0 pyfftw cython wget astroquery

export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
export FC=/opt/local/bin/gfortran

unset DYLD_LIBRARY_PATH
unset LIBRARY_PATH

cd /opt/anaconda3/envs/xcovspec/lib/python3.10/site-packages
wget https://www.astro.rug.nl/software/kapteyn/kapteyn-3.4.tar.gz
tar -xvf kapteyn-3.4.tar.gz
rm -f kapteyn-3.4.tar.gz
cd kapteyn-3.4
export CFLAGS="-Wno-error=int-conversion"
python setup.py install

cd ~/software/XcovSpec/
export PYTHONPATH=/Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}
