# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB3RD=$(PWD)/lib3rd/lib
INC3RD=$(PWD)/lib3rd/include

LIB_LAPACK_WRAPPER = liblapack_wrapper.a

CC    = gcc
CXX   = g++
F90   = gfortran
INC   =
LIBS  =
CLIBS = -lc++
DEFS  =

CXXFLAGS = -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -funroll-loops -fPIC
override LIBS += -L./lib -llapack_wrapper -L$(LIB3RD)
override INC  += -I./src -I$(INC3RD)

#
# select which version of BLAS/LAPACK use
#
USED_LIB=""
SYSTEMOPENBLAS=\/\/
ifeq ($(ATLAS),1)
  USED_LIB = LAPACK_WRAPPER_USE_ATLAS
endif
#
ifeq ($(MKL),1)
  USED_LIB = LAPACK_WRAPPER_USE_MKL
endif
#
ifeq ($(OPENBLAS),1)
  USED_LIB = LAPACK_WRAPPER_USE_OPENBLAS
endif
#
ifeq ($(LAPACK),1)
  USED_LIB = LAPACK_WRAPPER_USE_LAPACK
endif
#
ifeq ($(ACCELERATE),1)
  USED_LIB = LAPACK_WRAPPER_USE_ACCELERATE
endif
#
# if missig setup default
#
ifeq ($(USED_LIB), "")
ifeq (,$(wildcard .lapack_wrapper_config))
ifneq (,$(findstring Darwin, $(OS)))
  USED_LIB = LAPACK_WRAPPER_USE_ACCELERATE
else
  USED_LIB = LAPACK_WRAPPER_USE_LAPACK
endif
else
  USED_LIB = $(shell cat .lapack_wrapper_config)
endif
endif

$(shell echo "$(USED_LIB)" > .lapack_wrapper_config )
$(info $(USED_LIB))

#      # #    # #    # #    #
#      # ##   # #    #  #  #
#      # # #  # #    #   ##
#      # #  # # #    #   ##
#      # #   ## #    #  #  #
###### # #    #  ####  #    #

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN   = -Wall
  CC     = gcc
  CXX    = g++ -std=c++11 -pthread
  THREAD = LAPACK_WRAPPER_USE_THREAD
  #
  # activate C++11 for g++ >= 4.9
  #
  VERSION  = $(shell $(CC) -dumpversion)
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = ar rcs
  LIBSGCC = -lstdc++ -lm -pthread

ifneq (,$(findstring LAPACK_WRAPPER_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += -L$(LIB3RD) -Wl,-rpath,$(LIB3RD) -lopenblas -L$(FPATH)/../../.. -Wl,-rpath,$(LIB3RD)/../../..  -lgfortran
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ATLAS,$(USED_LIB)))
  ATLAS_PATH = /usr/lib/atlas-base
  ATLAS_LIBS = -llapack -lf77blas -lcblas -latlas -lgfortran
  override LIBS += -L$(ATLAS_PATH) $(ATLAS_LIBS) -Wl,-rpath,$(ATLAS_PATH)
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_MKL,$(USED_LIB)))
  # for MKL
  MKL_PATH = /opt/intel/mkl
  MKL_ARCH = intel64
  #MKL_LIB = -lmkl_tbb_thread -lmkl_rt -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_rt -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib/$(MKL_ARCH) -Wl,-rpath,$(MKL_PATH)/lib/$(MKL_ARCH) $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ACCELERATE,$(USED_LIB)))
  $(error error is "Accelerate is supported only on Darwin!")
endif

  #
  override INC  += -I/usr/include/eigen3 -I/usr/include/atlas/
endif

#######  #####  #     #
#     # #     #  #   #
#     # #         # #
#     #  #####     #
#     #       #   # #
#     # #     #  #   #
#######  #####  #     #

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
  CC       = clang
  CXX      = clang++
  VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
  #---------
  CXX     += -std=c++11 -stdlib=libc++
  THREAD   = LAPACK_WRAPPER_USE_THREAD
  #---------
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = libtool -static -o
  LIBSGCC = -lstdc++ -lm

ifneq (,$(findstring LAPACK_WRAPPER_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += -L$(LIB3RD) -Wl,-rpath,$(LIB3RD) -lopenblas -L$(FPATH)/../../.. -Wl,-rpath,$(LIB3RD)/../../..  -lgfortran
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ATLAS,$(USED_LIB)))
  $(error error is "Atlas is supported only on Linux!")
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_MKL,$(USED_LIB)))
  # for MKL
  MKL_PATH = /opt/intel/mkl
  #MKL_LIB = -lmkl_tbb_thread -lmkl_intel -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_intel -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib -Xlinker -rpath -Xlinker $(MKL_PATH)/lib $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ACCELERATE,$(USED_LIB)))
  override LIBS += -L./lib -llapack_wrapper -framework Accelerate
endif

  override INC += -I/usr/local/include/eigen3
endif

SRCS = \
src/lapack_wrapper/lapack_wrapper++.cc \
src/lapack_wrapper/lapack_wrapper.cc \
src/HSL/ma48_wrapper.cc \
src/HSL/ma57_wrapper.cc

OBJS = $(SRCS:.cc=.o)

SRCS_TESTS = \
src_tests/test1-small-factorization.cc \
src_tests/test2-Timing.cc \
src_tests/test3-BandedMatrix.cc \
src_tests/test4-BFGS.cc \
src_tests/test5-BLOCKTRID.cc \
src_tests/test6-EIGS.cc

OBJS_TESTS = $(SRCS_TESTS:.cc=.o)

#src/AlglinConfig.hh
DEPS = \
src/lapack_wrapper/TicToc.hh \
src/lapack_wrapper/lapack_wrapper++.hh \
src/lapack_wrapper/lapack_wrapper.hh \
src/lapack_wrapper/lapack_wrapper_config.hh \
src/lapack_wrapper/code/banded.hxx \
src/lapack_wrapper/code/blas.hxx \
src/lapack_wrapper/code/general.hxx \
src/lapack_wrapper/code/general_qr.hxx \
src/lapack_wrapper/code/sparse.cxx \
src/lapack_wrapper/code/sparse.hxx \
src/lapack_wrapper/code/symmetric.hxx \
src/lapack_wrapper/code/triangular.hxx \
src/lapack_wrapper/code/tridiagonal.hxx \
src/lapack_wrapper/code/wrapper.cxx \
src/lapack_wrapper/code/wrapper.hxx \
src/lapack_wrapper/HSL/hsl.h \
src/lapack_wrapper/HSL/hsl_solver.h \
src/lapack_wrapper/HSL/ma48_wrapper.h \
src/lapack_wrapper/HSL/ma57_wrapper.h

MKDIR = mkdir -p

.SUFFIXES:           # Delete the default suffixes
.SUFFIXES: .c .cc .o # Define our suffix list

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX = /usr/local

all: config lib $(OBJS_TESTS)
	mkdir -p bin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.o  $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test2-Timing              src_tests/test2-Timing.o               $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test3-BandedMatrix        src_tests/test3-BandedMatrix.o         $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test4-BFGS                src_tests/test4-BFGS.o                 $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test5-BLOCKTRID           src_tests/test5-BLOCKTRID.o            $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test6-EIGS                src_tests/test6-EIGS.o                 $(LIBS) $(LIBSGCC)

lib: config lib/$(LIB_LAPACK_WRAPPER)

.cc.o:
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

.c.o:
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/liblapack_wrapper.a: $(OBJS)
	$(MKDIR) lib
	$(AR) lib/liblapack_wrapper.a $(OBJS)

lib/liblapack_wrapper.dylib: $(OBJS)
	$(MKDIR) lib
	$(CXX) -shared -o lib/liblapack_wrapper.dylib $(OBJS)

lib/liblapack_wrapper.so: $(OBJS)
	$(MKDIR) lib
	$(CXX) -shared -o lib/liblapack_wrapper.so $(OBJS)

install_local:
	$(MKDIR) ./lib
	$(MKDIR) ./lib/include
	$(MKDIR) ./lib/include/lapack_wrapper
	$(MKDIR) ./lib/include/lapack_wrapper/code
	$(MKDIR) ./lib/include/lapack_wrapper/code++
	$(MKDIR) ./lib/include/HSL
	$(MKDIR) ./lib/include/sparse_tool
	$(MKDIR) ./lib/include/zstream
	cp -f -P src/lapack_wrapper/*.h*        ./lib/include/lapack_wrapper
	cp -f -P src/lapack_wrapper/code/*.h*   ./lib/include/lapack_wrapper/code
	cp -f -P src/lapack_wrapper/code++/*.h* ./lib/include/lapack_wrapper/code++
	cp -f -P src/HSL/*.h*                   ./lib/include/HSL
	cp -rf src/sparse_tool/*                ./lib/include/sparse_tool
	cp -rf src/zstream/*.h*                 ./lib/include/zstream
	cp -rf lib3rd/include/*                 ./lib/include

install: lib/$(LIB_LAPACK_WRAPPER)
	$(MKDIR) $(PREFIX)/include
	cp -rf lib/include/* $(PREFIX)/include
	cp -f -P lib/lib/*   $(PREFIX)/lib

config:
	rm -f src/lapack_wrapper/lapack_wrapper_config.hh
	sed 's/@@LAPACK_WRAPPER_USE@@/#define $(USED_LIB) 1/' <src/lapack_wrapper/lapack_wrapper_config.hh.tmpl | \
	sed 's/@@LAPACK_WRAPPER_THREAD@@/#define $(THREAD) 1/' | \
	sed 's/@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@/$(SYSTEMOPENBLAS)#define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1/' >src/lapack_wrapper/lapack_wrapper_config.hh

run:
	./bin/test1-small-factorization
	./bin/test2-Timing
	./bin/test3-BandedMatrix
	./bin/test4-BFGS
	./bin/test5-BLOCKTRID
	./bin/test6-EIGS

doc:
	doxygen

clean:
	rm -rf lib/liblapack_wrapper.* src/*.o src/*/*.o src/*/*/*.o src_tests/*.o
	rm -rf bin
