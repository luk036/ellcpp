# install dependencies
#Should match Config.py
echo Initializing...

#Create directories to store persistent data
# mkdir -p /workspace/conda
# mkdir -p /workspace/data

#Create a new env called arcw
# conda create --prefix /workspace/conda/arcw python=3.6 &&
# echo "conda activate /workspace/conda/arcw" >> ~/.bashrc &&
export CONDA_PREFIX=/opt/conda
export PATH=$CONDA_PREFIX/bin:$PATH &&
source ~/.bashrc
export SHELL=/bin/bash

#Install conda packages for to run jupyterlab
# hash -r
# conda config --set always_yes yes --set changeps1 no
# conda update -q conda
# conda info -a
# conda install -y ninja
# conda install -y -c conda-forge catch2
# conda install -y -c conda-forge libboost
# conda install -y -c conda-forge cppcheck

# conda install -y -c conda-forge jupyterlab
# conda install -y -c conda-forge beakerx
# conda install -y -c conda-forge xeus-cling
# conda install -y -c conda-forge xeus-python

# Some extra packages for your environments
# conda install -y -c anaconda pyodbc
# conda install -y -c conda-forge xeus-cling
# conda install -y -c conda-forge xeus-python
# conda install -y -c conda-forge python-graphviz
# conda install -y -c conda-forge tensorflow
# conda install -y -c conda-forge keras
# conda install -y -c conda-forge xtl
# conda install -y -c conda-forge openblas
# conda install -y -c conda-forge gdal
# conda install -y -c conda-forge util-linux
# conda install -y -c conda-forge libtiff
# conda install -y -c conda-forge libgdal
# conda install -y -c pytorch pytorch
# conda install -y -c conda-forge dask
# conda install -y -c conda-forge dash
# conda install -y -c conda-forge dash-table
# conda install -y -c conda-forge rx
# conda install -y -c conda-forge dash-core-components
# conda install -y -c conda-forge cassandra-driver
# conda install -y -c ranaroussi yfinance

#Notes
#To run jupyter lab : jupyter lab --ip=0.0.0.0 --allow-root
echo Done...