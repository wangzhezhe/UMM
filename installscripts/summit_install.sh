#!/bin/bash

set -e

# refer to https://github.com/CODARcode/MGARD/blob/master/build_scripts/build_mgard_cuda_summit.sh
# for installing the mgard
# put the gpu installation as needed

module load gcc
module load cmake

HERE=`pwd`
build_jobs=8
source ./settings.sh
SOFTWARE_SRC_DIR="$HERE/src"
SOFTWARE_BUILD_DIR="$HERE/build"
SOFTWARE_INSTALL_DIR="$HERE/install"

if [ ! -d $SOFTWARE_SRC_DIR ];
    then
    mkdir $SOFTWARE_SRC_DIR
fi
echo "SOFTWARE_SRC_DIR: $SOFTWARE_SRC_DIR"


echo "====> Installing zstd"
ZSTD_SRC_DIR="$SOFTWARE_SRC_DIR/zstd"
ZSTD_BUILD_DIR="$SOFTWARE_BUILD_DIR/zstd"
ZSTD_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/zstd"

if [ -d $ZSTD_INSTALL_DIR ]; then

echo "====> skip, $ZSTD_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else

    if [ ! -d $ZSTD_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone $ZSTD_REPO
    cd $ZSTD_SRC_DIR
    fi


echo "**** Configuring ZSTD"
cmake -S ${ZSTD_SRC_DIR}/build/cmake -B ${ZSTD_BUILD_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=${ZSTD_INSTALL_DIR} \
  -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc

echo "**** Building ZSTD"
cmake --build ${ZSTD_BUILD_DIR} -j${build_jobs}
echo "**** Installing ZSTD"
cmake --install ${ZSTD_BUILD_DIR}
fi



echo "====> Installing protobuf"
PROTO_SRC_DIR="$SOFTWARE_SRC_DIR/protobuf"
PROTO_BUILD_DIR="$SOFTWARE_BUILD_DIR/protobuf"
PROTO_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/protobuf"

if [ -d $PROTO_INSTALL_DIR ]; then

echo "====> skip, $PROTO_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else

    if [ ! -d $PROTO_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone -b v3.19.4 $PROTOBUF_REPO
    cd $PROTO_SRC_DIR
    git submodule update --init --recursive
    fi


echo "**** Configuring protobuf"
cmake -S ${PROTO_SRC_DIR}/cmake -B ${PROTO_BUILD_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -Dprotobuf_BUILD_TESTS=OFF \
  -Dprotobuf_BUILD_SHARED_LIBS=ON \
  -DCMAKE_INSTALL_PREFIX=${PROTO_INSTALL_DIR} \
  -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc

echo "**** Building ZSTD"
cmake --build ${PROTO_BUILD_DIR} -j${build_jobs}
echo "**** Installing ZSTD"
cmake --install ${PROTO_BUILD_DIR}
fi


echo "====> Installing mgard"
MGARD_SRC_DIR="$SOFTWARE_SRC_DIR/MGARD"
MGARD_BUILD_DIR="$SOFTWARE_BUILD_DIR/MGARD"
MGARD_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/MGARD"

if [ -d $MGARD_INSTALL_DIR ]; then

echo "====> skip, $MGARD_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else

    if [ ! -d $MGARD_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone $MGARD_REPO
    cd $MGARD_SRC_DIR
    fi

echo "**** Configuring MGARD"
cmake -S ${MGARD_SRC_DIR} -B ${MGARD_BUILD_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=${MGARD_INSTALL_DIR} \
  -DCMAKE_PREFIX_PATH="${ZSTD_INSTALL_DIR}/lib/cmake/zstd;${PROTO_INSTALL_DIR}" \
  -Dzstd_DIR=${ZSTD_INSTALL_DIR}/lib64/cmake/zstd \
  -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc

echo "**** Building MGARD"
cmake --build ${MGARD_BUILD_DIR} -j${build_jobs}
echo "**** Installing MGARD"
cmake --install ${MGARD_BUILD_DIR}
fi

echo "====> Installing VTK"
VTK_SRC_DIR="$SOFTWARE_SRC_DIR/vtk"
VTK_BUILD_DIR="$SOFTWARE_BUILD_DIR/vtk"
VTK_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/vtk"

if [ -d $VTK_INSTALL_DIR ]; then

echo "====> skip, $VTK_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else

    if [ ! -d $VTK_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone $VTK_REPO
    cd $VTK_SRC_DIR
    git submodule update --init --recursive 
    fi

echo "**** Configuring VTK"
cmake -S ${VTK_SRC_DIR} -B ${VTK_BUILD_DIR} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=${VTK_INSTALL_DIR} \
  -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc

echo "**** Building VTK"
cmake --build ${VTK_BUILD_DIR} -j${build_jobs}
echo "**** Installing VTK"
cmake --install ${VTK_BUILD_DIR}
fi


echo "====> Building UMM"
UMM_SRC_DIR=$HERE/../
# use the install dir as the build dir
UMM_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/UMM"

if [ -d $UMM_INSTALL_DIR ]; then
    echo "====> skip, $UMM_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else

    cmake -B ${UMM_INSTALL_DIR} -S ${UMM_SRC_DIR} \
    -DVTK_DIR=${VTK_INSTALL_DIR}/lib64/cmake/vtk-9.2 \
    -DCMAKE_PREFIX_PATH="${ZSTD_INSTALL_DIR}/lib/cmake/zstd;${PROTO_INSTALL_DIR}" \
    -Dmgard_DIR=${MGARD_INSTALL_DIR}/lib64/cmake/mgard \
    -Dzstd_DIR=${ZSTD_INSTALL_DIR}/lib64/cmake/zstd \
    -DBUILD_TESTING=OFF \
    -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
    
    cd $HERE

    # build and install
    echo "**** Building UMM"
    cmake --build ${UMM_INSTALL_DIR} -j${build_jobs}
fi