OpenBlas 0.3.17 compiled for windows 64 bit

compiled using MSYS2

pacman -S mingw-w64-x86_64-toolchain
pacman -S mingw-w64-x86_64-openmp

export COMMON_OPT="-O3 -ftree-vectorize -fprefetch-loop-arrays --param prefetch-latency=300"
export CFLAGS="-O3 -ftree-vectorize -fprefetch-loop-arrays --param prefetch-latency=300"
export FCOMMON_OPT="-O3 -ftree-vectorize -fprefetch-loop-arrays --param prefetch-latency=300"
export FCFLAGS="-O3 -ftree-vectorize -fprefetch-loop-arrays --param prefetch-latency=300"

NO_LAPACKE=1

# compile threaded version with OPENMP
make -j DYNAMIC_ARCH=1 USE_OPENMP=1 NUM_PARALLEL=16 GEMM_MULTITHREAD_THRESHOLD=50 BUILD_RELAPACK=1 CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran AR=x86_64-w64-mingw32-gcc-ar RANLIB=x86_64-w64-mingw32-gcc-ranlib

# compile threaded version with THREAD
make -j DYNAMIC_ARCH=1 GEMM_MULTITHREAD_THRESHOLD=50 BUILD_RELAPACK=1 CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran AR=x86_64-w64-mingw32-gcc-ar RANLIB=x86_64-w64-mingw32-gcc-ranlib

# compile non threaded version
make -j DYNAMIC_ARCH=1 CC=x86_64-w64-mingw32-gcc FC=x86_64-w64-mingw32-gfortran AR=x86_64-w64-mingw32-gcc-ar RANLIB=x86_64-w64-mingw32-gcc-ranlib BINARY=64 INTERFACE=64 NO_AFFINITY=1 NO_WARMUP=1 USE_OPENMP=0 USE_THREAD=0 USE_LOCKING=1

# compile native via conda+FLANG
https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio

conda update -n base conda
conda config --add channels conda-forge
conda install -y cmake flang clangdev perl libflang ninja

"c:/Program Files (x86)/Microsoft Visual Studio/2017/Community/VC/Auxiliary/Build/vcvars64.bat"

set "LIB=%CONDA_PREFIX%\Library\lib;%LIB%"
set "CPATH=%CONDA_PREFIX%\Library\include;%CPATH%"
mkdir build
cd build

cmake .. -G "Ninja" -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_C_COMPILER=clang-cl -DCMAKE_Fortran_COMPILER=flang -DCMAKE_MT=mt -DBUILD_WITHOUT_LAPACK=no -DNOFORTRAN=0 -DDYNAMIC_ARCH=ON -DCMAKE_BUILD_TYPE=Release

cmake --build . --config Release

cmake --install . --prefix c:\opt -v
