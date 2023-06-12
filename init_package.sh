#!/bin/bash


#Setting up environment

# Require sudo privilegies
if [ $UID == 0 ]; then
    echo "Please run this script without sudo:"
    echo "$0 $*"
    exit 1
fi

make_install_build_of_type(){
    echo
    echo
    echo
    echo
    echo "Cmake and install build of type " $1
    echo
    echo
    echo
    echo
    mkdir -p $1 && cd $1
    cmake ../.. -DCMAKE_BUILD_TYPE=$1
    make
    sudo make install
}

echo
echo
echo
echo
echo "setting up build system"
echo
echo
echo
echo
mkdir -p build && cd build

make_install_build_of_type Debug
cd ..
make_install_build_of_type Release


#make doc

