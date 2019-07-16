# VFFTPACK-CMake (https://github.com/dmey/vfftpack-cmake).
# Copyright (C) 2019 D. Meyer.

cmake_minimum_required(VERSION 3.0)
project(VFFTPACK Fortran)

include(CTest)

option(USE_REAL8 "Use real8 instead of real4" ON)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -fbacktrace  -Wextra -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all")
    if (USE_REAL8)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
    endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check bounds -traceback -fpe-all0 -check all")
    if (USE_REAL8)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
    endif()
    else()
    message(FATAL_ERROR "Only gfortran and ifort compilers are supported")
endif()

file(GLOB_RECURSE vfft_src "src/*.f")

# Use "object libraries" to avoid the duplicated compilation
# of the object files (https://stackoverflow.com/a/29824424).
# This is the "object library" target: compiles the sources only once
add_library(objlib OBJECT ${vfft_src})

# Shared libraries need PIC
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# Shared and static libraries built from the same object files
add_library(vfftpack_shared SHARED $<TARGET_OBJECTS:objlib>)
add_library(vfftpack_static STATIC $<TARGET_OBJECTS:objlib>)