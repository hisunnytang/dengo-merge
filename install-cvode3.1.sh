set -ex
tar -xvzf cvode-3.1.0.tar.gz
cd cvode-3.1.0/builddir
cmake -DCMAKE_INSTALL_PREFIX=~/virtualenv/cvode_instdir -DKLU_INCLUDE_DIR=~/virtualenv/suitesparse/include -DKLU_LIBRARY_DIR=~/virtualenv/suitesparse/KLU/lib -DCMAKE_C_FLAGS=-fPIC -DKLU_ENABLE=ON ../
make && make install
