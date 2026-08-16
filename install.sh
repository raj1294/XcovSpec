conda create --name XcovSpec python==3.10 -c https://cxc.cfa.harvard.edu/conda/ciao ciao caldb_main
source activate XcovSpec
export PYTHONPATH=/opt/homebrew/lib/python3.10/site-packages:${PYTHONPATH}
python -m pip install astropy scipy scikit-learn stingray numpy==2.0.0 numba==0.67.0 cython 

