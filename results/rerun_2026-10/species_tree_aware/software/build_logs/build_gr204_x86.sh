#!/bin/bash
B=/private/tmp/claude-501/-Users-spencergray-Desktop-pif1-foldtree/c1b55069-62c6-4237-9ba9-4a315cadb526/scratchpad/build
source ~/miniforge3/etc/profile.d/conda.sh; conda activate strecon_x86
export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.sdk CONDA_BUILD_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.sdk
cd $B/generax204 && rm -rf build_x86 && mkdir build_x86 && cd build_x86
arch -x86_64 cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_OSX_ARCHITECTURES=x86_64 -DCMAKE_OSX_SYSROOT=$SDKROOT -DBISON_EXECUTABLE=$(which bison) -DFLEX_EXECUTABLE=$(which flex) > $B/gr204_x86_cmake.log 2>&1 || { echo cmake failed; exit 1; }
arch -x86_64 make -j4 generax > $B/gr204_x86_make.log 2>&1; echo "make exit $?"
file bin/generax
