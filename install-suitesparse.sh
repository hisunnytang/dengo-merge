set -ex
git clone https://github.com/jluttine/suitesparse.git
cd suitesparse
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS && make && make install PREFIX="/home/travis/build/hisunnytang/dengo-merge/suitesparse/OpenBLAS/Install"
cd ../
export LAPACK=""
export BLAS="-L/home/travis/build/hisunnytang/dengo-merge/suitesparse/OpenBLAS/Install/lib -lopenblas"
cd SuiteSparse_config && make && cd ../
cd AMD && make && cd ../
cd CAMD && make && cd ../
cd CCOLAMD && make && ../
cd metis-5.1.0 && make config && make && make install PREFIX="/home/travis/build/hisunnytang/dengo-merge/suitesparse/OpenBLAS/Install/lib" && cd ../
cd COLAMD && make && cd ../
cd BTF && make && cd ../
cd KLU && make
