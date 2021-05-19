git clone https://github.com/DrTimothyAldenDavis/SuiteSparse
cd SuiteSparse
make purge
# make distclean
export CUDA=no
export BLAS="-L/home/cwagner4/miniconda/lib -lblas"
export LAPACK="-L/home/cwagner4/miniconda/lib -llapack"

cd SuiteSparse_config && make CUDA=no &> /dev/null && cd ../
make metis CUDA=no &> /dev/null
cd AMD && make CUDA=no &> /dev/null && cd ../
cd CAMD && make CUDA=no &> /dev/null && cd ../
cd CCOLAMD && make CUDA=no &> /dev/null && cd ../
cd COLAMD && make CUDA=no &> /dev/null && cd ../
cd BTF && make CUDA=no &> /dev/null && cd ../
cd CHOLMOD && make CUDA=no && cd ../
cd KLU && make CUDA=no
