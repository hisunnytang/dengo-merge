wget https://download.chiantidatabase.org/CHIANTI_10.0.1_database.tar.gz
mkdir CHIANTI_Data/
tar -xvzf CHIANTI_10.0.1_database.tar.gz -C CHIANTI_data
export XUVTOP=$(TRAVIS_BUILD_DIR)/CHIANTI_data
pip3 install ChiantiPy

