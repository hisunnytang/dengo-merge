set -ex
tar -xvzf cvode-3.1.0.tar.gz
cd cvode-3.1.0 && rm -rf builddir && mkdir builddir
cd builddir
cvode_inst="${TRAVIS_BUILD_DIR}/cvode_instdir"
suitesparse_inst="${TRAVIS_BUILD_DIR}/suitesparse"
cmake -DCMAKE_INSTALL_PREFIX=$cvode_inst -DKLU_INCLUDE_DIR="${suitesparse_inst}/include" -DKLU_LIBRARY_DIR="${suitesparse}/KLU/lib" -DCMAKE_C_FLAGS=-fPIC -DKLU_ENABLE=ON ../
make && make install
