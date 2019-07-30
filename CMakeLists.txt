# Disable in-source builds to prevent source tree corruption.
IF (" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
	message(FATAL_ERROR "FATAL: In-source builds are not allowed.
	You should create a separate directory for build files.")
ENDIF()

# Set minimum cmake version + cmake settings
cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 OLD)

message(STATUS "Configure F5c")

# Set project name
project (F5C LANGUAGES CXX C)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -rdynamic -Wall -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -rdynamic -Wall -O2 -std=c++11 ")

# CMAKE_EXE_LINKER_FLAGS
# CMAKE_MODULE_LINKER_FLAGS
# CMAKE_SHARED_LINKER_FLAGS
# CMAKE_STATIC_LINKER_FLAGS

add_subdirectory(${CMAKE_SOURCE_DIR}/src)

