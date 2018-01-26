# assumes conda
conda create -n glicko --file requirements.txt -c conda-forge -y
source activate glicko
cd glicko2 && python setup.py install
