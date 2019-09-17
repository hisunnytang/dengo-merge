set -ex
git clone https://github.com/jluttine/suitesparse.git
cd suitesparse
git clone https://github.com/xianyi/OpenBLAS.git
openblas_inst="${TRAVIS_BUILD_DIR}/suitesparse/OpenBLAS/Install"
cd OpenBLAS && make && make install PREFIX=$openblas_inst
cd ../
export LAPACK=""
export BLAS="-L${openblas_inst}/lib -lopenblas"
cd SuiteSparse_config && make && cd ../
cd AMD && make && cd ../
cd CAMD && make && cd ../
cd CCOLAMD && make && cd ../
suitesparse_inst="${TRAVIS_BUILD_DIR}/suitesparse"
cd metis-5.1.0 && make config prefix=$suitesparse && make && make install && cd ../
cd COLAMD && make && cd ../
cd BTF && make && cd ../
cd KLU && make
