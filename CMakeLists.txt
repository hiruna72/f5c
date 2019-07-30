# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.")
endif()

# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 OLD)

# Specify the project, and set the use of c++ and c
project(libtwentyface LANGUAGES CXX C)

set(HTS_VERSION 1.9)
set(HDF5_VERSION 1.10.4)
set(HDF5_MAJOR_MINOR 1.10)

# configure step (the first step)
execute_process(
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND curl -o ./htslib.tar.bz2 -L https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2
    COMMAND tar -xf htslib.tar.bz2  
)

execute_process(
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND curl -o ./hdf5.tar.bz2 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_MAJOR_MINOR}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.bz2
    COMMAND tar -xf hdf5.tar.bz2  
)

# second step (build)
add_custom_target(htslib ALL
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/htslib-${HTS_VERSION}
    COMMAND ./configure --prefix=`pwd`/../ --enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no
    COMMAND make -j4
    COMMAND make install
)

add_custom_target(hdf5 ALL
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/hdf5-${HDF5_VERSION}
    COMMAND ./configure --prefix=`pwd`/../
    COMMAND make -j4
    COMMAND make install
)