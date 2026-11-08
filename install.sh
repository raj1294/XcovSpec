conda create --name xcovspec python==3.10
conda install --name xcovspec -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge ciao caldb_main ciao-contrib sherpa llvmlite mpmath ds9
source activate xcovspec
export CC=/opt/local/bin/gcc
export CXX=/opt/local/bin/g++
export FC=/opt/local/bin/gfortran
export PYTHON=/opt/local/bin/python
export PATH=/opt/anaconda3/bin:/opt/local/bin:${PATH}
python -m pip install astropy scipy scikit-learn stingray numpy==1.26.0 numba cython wget astroquery
cd /Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages
wget https://www.astro.rug.nl/software/kapteyn/kapteyn-3.4.tar.gz
tar -xvf kapteyn-3.4.tar.gz
rm -f kapteyn-3.4.tar.gz
cd kapteyn-3.4
export CFLAGS="-Wno-error=int-conversion"
python setup.py install
cd ~/software/XcovSpec/
python -m pip install numpy==2.0.0
export PYTHONPATH=/Users/raj1294/anaconda3/envs/xcovspec/lib/python3.10/site-packages:${PYTHONPATH}
