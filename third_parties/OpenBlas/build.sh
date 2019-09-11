export OPENBLAS_NUM_THREADS=0
export GOTO_NUM_THREADS=0
export OMP_NUM_THREADS=0
make DYNAMIC_ARCH=1 BINARY=64 BINARY64=1

#USE_THREAD
#USE_OPENMP
make DYNAMIC_ARCH=1 BINARY=64 BINARY64=1 \
     HOSTCC=gcc CC=x86_64-w64-mingw32-gcc \
     FC=x86_64-w64-mingw32-gfortran \
     CFLAGS="-static-libgcc -static-libstdc++ -static -g0" \
     FFLAGS="-static"

make clean

make DEBUG=1 DYNAMIC_ARCH=1 BINARY=64 BINARY64=1 \
     HOSTCC=gcc CC=x86_64-w64-mingw32-gcc \
     FC=x86_64-w64-mingw32-gfortran \
     CFLAGS="-static-libgcc -static-libstdc++ -static -ggdb" \
     FFLAGS="-static"
