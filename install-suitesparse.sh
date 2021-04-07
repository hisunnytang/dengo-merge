set -ex
git clone https://github.com/DrTimothyAldenDavis/SuiteSparse
cd SuiteSparse
make purge
make distclean
export LAPACK=""
export BLAS=""
cd SuiteSparse_config && make CUDA=no && cd ../
make metis CUDA=no
cd AMD && make CUDA=no && cd ../
cd CAMD && make CUDA=no && cd ../
cd CCOLAMD && make CUDA=no && cd ../
cd COLAMD && make CUDA=no && cd ../
cd BTF && make CUDA=no && cd ../
cd KLU && make CUDA=no
