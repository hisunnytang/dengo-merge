set -ex
git clone https://github.com/jluttine/suitesparse.git
cd suitesparse/SuiteSparse_config && make && cd ../
cd AMD && make && cd ../
cd BTF && make && cd ../
cd KLU && make
