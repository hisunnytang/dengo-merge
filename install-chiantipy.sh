#wget https://download.chiantidatabase.org/CHIANTI_10.0.1_database.tar.gz
mkdir CHIANTI_Data
tar -xvzf CHIANTI_10.0.1_database.tar.gz -C CHIANTI_Data
export XUVTOP=$TRAVIS_BUILD_DIR/CHIANTI_Data
pip install ChiantiPy

