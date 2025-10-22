#!/bin/bash
#
eval "$(conda shell.bash hook)"
#
conda create -n genomeclass -y
conda activate genomeclass
conda install -c conda-forge python=3.13.0 -y
conda install -c conda-forge ncbi-datasets-cli -y
conda install -c bioconda geco3 -y
conda install -c bioconda jarvis3 -y
conda install -c bioconda spades -y
#
pip install scikit-learn
pip install nltk
pip install matplotlib
pip install pandas
pip install seaborn
pip install numpy
pip install xgboost
pip install shap
pip install ipython
#
conda activate base
#
