#!/bin/bash 

# Defining path variables
DENGO_INSTALL_PATH=$PWD
CVODE_PATH=$DENGO_INSTALL_PATH/cvode_instdir
HDF5_PATH=$DENGO_INSTALL_PATH/hdf5_install
SUITESPARSE_PATH=$DENGO_INSTALL_PATH/SuiteSparse
HDF5_DIR=$DENGO_INSTALL_PATH/hdf5_install

# Exporting path variables
echo "export DENGO_INSTALL_PATH=$PWD" >> ~/.bashrc
echo "export TRAVIS_BUILD_DIR=$DENGO_INSTALL_PATH" >> ~/.bashrc #To remain consistent with naming conventions in Travis
echo "export CVODE_PATH=$DENGO_INSTALL_PATH/cvode_instdir" >> ~/.bashrc
echo "export HDF5_PATH=$DENGO_INSTALL_PATH/hdf5_install" >> ~/.bashrc
echo "export SUITESPARSE_PATH=$DENGO_INSTALL_PATH/SuiteSparse" >> ~/.bashrc
echo "export HDF5_DIR=$DENGO_INSTALL_PATH/hdf5_install" >> ~/.bashrc

# Checking if the user would like to provide a different HDF5 path
echo "HDF5_DIR=$HDF5_DIR"
while true; do
    read -p "Would you like to specify a different HDF5 path? " yn
    case $yn in
        [Yy]* ) read -p "HDF5_PATH=" path; read -p "New HDF5_DIR=$path. Is this correct? " yn
            case $yn in
                [Yy]* ) echo "export HDF5_PATH=$path"; break;;
                [Nn]* ) ;;
                * ) echo "Please answer yes or no" ;;
            esac ;;
        [Nn]* ) path=$HDF5_DIR; break ;;
        * ) echo "Please answer yes or no" ;;
    esac
done

# Appending variables to library path
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DENGO_INSTALL_PATH/lib:$HDF5_PATH/lib:$SUITESPARSE_PATH/lib:$CVODE_PATH/lib" >> ~/.bashrc

source ~/.bashrc

# Running separate install scripts
./install-suitesparse.sh &> suitesparse_install.out
./install-cvode.sh &> cvode_install.out
./install-chiantipy.sh &> chiantipy_install.out
./install-yt.sh &> yt_install.out

if [ "$HDF5_DIR" = "$path" ]; then
    ./install-hdf5.sh &> hdf5_install.out 
fi

# Install Python packages
pip install -r requirements.txt
pip install h5py
pip install sympy
pip install yt
pip install -e .