# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

CC    = gcc
CXX   = g++
F90   = gfortran
INC   =
LIBS  =
CLIBS = -lc++
DEFS  =

CXXFLAGS = -O2 -funroll-loops -fPIC
override INC  += -I./src -Ilib3rd/include
override LIBS += -Llib3rd/lib -Llib3rd/dll

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
  USED_LIB = LAPACK_WRAPPER_USE_OPENBLAS
endif
else
  USED_LIB = $(shell cat .lapack_wrapper_config)
endif
endif

$(shell echo "$(USED_LIB)" > .lapack_wrapper_config )
$(info $(USED_LIB))

SRCS = \
src/lapack_wrapper/lapack_wrapper++.cc \
src/lapack_wrapper/lapack_wrapper.cc \
src/HSL/hsl_solver.cc \
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

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
include Makefile_linux.mk
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
include Makefile_mingw.mk
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
include Makefile_osx.mk
endif

MKDIR = mkdir -p

.SUFFIXES:           # Delete the default suffixes
.SUFFIXES: .c .cc .o # Define our suffix list

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX = /usr/local

tests: all_libs $(OBJS_TESTS)
	mkdir -p bin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.o  $(ALL_LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test2-Timing              src_tests/test2-Timing.o               $(ALL_LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test3-BandedMatrix        src_tests/test3-BandedMatrix.o         $(ALL_LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test4-BFGS                src_tests/test4-BFGS.o                 $(ALL_LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test5-BLOCKTRID           src_tests/test5-BLOCKTRID.o            $(ALL_LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test6-EIGS                src_tests/test6-EIGS.o                 $(ALL_LIBS) $(LIBSGCC)

.cc.o:
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

.c.o:
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

install_local: all_libs
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

install: all_libs
	$(MKDIR) $(PREFIX)/include
	cp -rf lib/include/* $(PREFIX)/include
	cp -f -P lib/lib/*   $(PREFIX)/lib

config:
	rm -f src/lapack_wrapper/lapack_wrapper_config.hh
	sed 's/@@LAPACK_WRAPPER_USE@@/#define $(USED_LIB) 1/' <src/lapack_wrapper/lapack_wrapper_config.hh.tmpl | \
	sed 's/@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@/$(SYSTEMOPENBLAS)#define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1/' >src/lapack_wrapper/lapack_wrapper_config.hh

run: tests
	./bin/test1-small-factorization
	./bin/test2-Timing
	./bin/test3-BandedMatrix
	./bin/test4-BFGS
	./bin/test5-BLOCKTRID
	./bin/test6-EIGS

doc:
	doxygen

clean:
	rm -rf lib/* src/*.o src/*/*.o src/*/*/*.o src_tests/*.o src/*.obj src/*/*.obj src/*/*/*.obj src_tests/*.obj
	rm -rf bin
