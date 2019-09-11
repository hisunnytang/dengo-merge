set -ex
git clone https://github.com/jluttine/suitesparse.git
cd suitesparse/SuiteSparse_config && make && cd ../
cd suitesparse/AMD && make && cd ../
cd suitesparse/COLAMD && make && cd ../
cd suitesparse/BTF && make && cd ../
cd suitesparse/KLU && make
