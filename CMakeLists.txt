# Disable in-source builds to prevent source tree corruption.
if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.")
endif()

# choose the architecture 
set(DEPLOY_PLATFORM "unknown" CACHE STRING "choose the instruction set architecture")
set_property(CACHE DEPLOY_PLATFORM PROPERTY STRINGS unknown x86 arm64-v8a armeabi-v7a)
IF (${DEPLOY_PLATFORM} STREQUAL "unknown")
	message(FATAL_ERROR "choose the DEPLOY_PLATFORM") 
	return() # This is to stop proceeding further and to stop opencv getting set to the default ANDROID_ABI
ENDIF()


# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 OLD)

# Specify the project, and set the use of c++ and c
project(F5C LANGUAGES CXX C)

set(HTS_VERSION 1.9)
set(HDF5_VERSION 1.10.4)
set(HDF5_MAJOR_MINOR 1.10)

message(STATUS "Configure F5c")

# Set project name
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fPIC -rdynamic -Wall -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -rdynamic -Wall -O2 -std=c++11 ")


add_subdirectory(${CMAKE_SOURCE_DIR}/src)
# add_subdirectory(${CMAKE_SOURCE_DIR}/examples)
