# math_tools

## Overview
This package contains the mathematical functions related to the Lie algebra and the Chebyshev differentiation.
Specifically, it implements all the transformation and operation on the SO(3) and SE(3) groups.
It also contains everything related to the Chebyshev points such as the grid and the differentiation matrix.
It provides a generic library that allows to integrate linear non-homogeneous ordinary differential equations.

The details about notations and definition are available in the compiled PDF.

If you are reading from GitHub <a href="Latex/math_tools.pdf">Click Here</a>.

If you are reading from Doxygen <a href="../../../Latex/math_tools.pdf">Click Here</a>


# Build and Install

##  Requirements

This package REQUIRES to have the mathematical library Eigen and Boost installed.

For Eigen, we need the version 3.4 or newest version.The latest library release can be found at the <a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">Boost Main Page</a>.

For Boost, the install procedure can be found at the <a href="https://www.boost.org/doc/libs/1_80_0/more/getting_started/unix-variants.html">Boost Main Page</a>.

## Suggested external libraries

In order to run tests and benchmarks this package uses Google Test and Google Benchmak.

The installation procedure for Google Test can be found at <a href="https://github.com/google/googletest.git">Google Test Main Page</a>.

The installation procedure for Google Benchmak can be found at <a href="https://github.com/google/benchmark.git">Google Benchmak Main Page</a>.

In order to generate documentation we use Doxygen; Doxygen can be installed following the procedure at the <a href="https://www.doxygen.nl/manual/install.html">Doxygen Installation Page</a>

## Installation procedure

Once the needed library is installed, you can download this repository in the preferred way.

Once inside the downloaded folder it is suffient to run the script init_package.sh running the following in the project folder.

    ./init_package.sh


### Build and install manually

If something goes wrong with the automatic installation, you can follow this procedure.
First, open a terminal in the package folder (same level as the README.md)

    mkdir -p build && cd build

    mkdir -p Debug && cd Debug

You can now build the code in Debug mode

    cmake ../.. -DCMAKE_BUILD_TYPE=Debug
    make

Intallation will be executed by

    sudo make install

You now step out of the Debug folder

    cd ..

You can now build the code in Release mode

    cmake ../.. -DCMAKE_BUILD_TYPE=Release
    make

Intallation will be executed by

    sudo make install

You can now build the documentation

    make doc

