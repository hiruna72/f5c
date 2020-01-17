#!/bin/sh

# for armeabi-v7a and arm64-v8a cross compiling
toolchain_file="set path to /android-sdk-linux/ndk-bundle/build/cmake/android.toolchain.cmake"

# terminate script
die(){
	echo "$1"
	exit 1
}


cd src
touch config.h
echo "#define HAVE_HDF5_H 1" > config.h

# to create a f5c library
sed -i 's/int main(int argc/int init_f5c(int argc/g' main.c || die "main not found"
touch interface.h
echo "int init_f5c(int argc, char *argv[]);" > interface.h

cd ..
mkdir -p build
rm -rf build
mkdir build
cd build

# for architecture x86
 # cmake .. -DDEPLOY_PLATFORM=x86
 # make -j 8

# # for architecture armeabi-V7a
# cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-21 -DDEPLOY_PLATFORM:STRING="armeabi-v7a" -DANDROID_ABI="armeabi-v7a"
# ninja

# # for architecture arm66-v8a
cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-23 -DDEPLOY_PLATFORM:STRING="arm64-v8a" -DANDROID_ABI="arm64-v8a"
ninja
