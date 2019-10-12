wget http://www.chiantidatabase.org/download/CHIANTI_9.0.1_data.tar.gz
tar -xvzf CHIANTI_9.0.1_data.tar.gz -C CHIANTI_data
export XUVTOP=$(TRAVIS_BUILD_DIR)/CHIANTI_data
pip install ChiantiPy

