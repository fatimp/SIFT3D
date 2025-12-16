# SIFT3D

This is a refactored version of SIFT3D by Blaine Rister et al. Changes include a
refactored and more programmer-friendly API, some bug fixes and code
cleanup. The original readme is in README-OLD.md file.

[Documentation](https://fatimp.github.io/SIFT3D)

## Dependencies

* LAPACK
* cmake (for building)
* nifticlib (for examples, optional)
* doxygen (for documentation, optional)

## Building instructions:

~~~~
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DWITH_EXAMPLES=ON ( or =OFF ) ..
make
make doc
make install
~~~~
