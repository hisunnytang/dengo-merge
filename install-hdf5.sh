wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
tar -xvzf hdf5-1.10.5.tar.gz
cd hdf5-1.10.5
./configure --prefix="${TRAVIS_BUILD_DIR}/hdf5_install"
make &> hdf5build.out
make install %> hdf5build.out
