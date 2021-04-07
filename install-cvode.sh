set -ex
# wget https://github.com/LLNL/sundials/releases/download/v5.7.0/cvode-5.7.0.tar.gz
# tar -xvzf cvode-5.7.0.tar.gz
cd cvode-5.7.0 && rm -rf builddir && mkdir builddir
cd builddir
cvode_inst="${TRAVIS_BUILD_DIR}/cvode_instdir"
suitesparse_inst="${TRAVIS_BUILD_DIR}/SuiteSparse"
cmake -DCMAKE_INSTALL_PREFIX=$cvode_inst -DKLU_INCLUDE_DIR="${suitesparse_inst}/include" -DKLU_LIBRARY_DIR="${suitesparse_inst}/lib" -DCMAKE_C_FLAGS=-fPIC -DENABLE_KLU=ON ../
make && make install
