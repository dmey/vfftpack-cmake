# VFFTPACK-CMake

This repository is a fork of [VFFTPACK Fortran 77 library version 2.1](http://www.netlib.no/netlib/vfftpack) with added CMake support.

For more information about the original library, see [doc](doc).


## Build instructions

### Prerequisites
The following software are required to build VFFTPACK-CMake from source:

- Linux or macOS
- Intel or GNU Fortran compiler
- [CMake](https://git-scm.com/) version 3.0, or above
- [Git](https://git-scm.com/) (_Optional_ if cloning)


### Download

Clone or download VFFTPACK-CMake from the GitHub repository at https://github.com/dmey/vfftpack-cmake.

### Build and tests

Go into `vfftpack-cmake` and run the following commands from your command-line interface:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_REAL8=ON ..
make # Creates shared and dynamic library
```

## Support

Please note that VFFTPACK-CMake is not being actively developed or supported.