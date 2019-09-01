WARN   = -Wall
CC     = gcc
CXX    = g++ -std=c++11 -pthread
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
  override INC  += -I/usr/include/atlas/
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

ALL_LIBS = $(LIBS) -Llib -llapack_wrapper_linux

all_libs: config lib/liblapack_wrapper_linux.so lib/liblapack_wrapper_linux_static.a

lib/liblapack_wrapper_linux.so: $(OBJS)
	@$(MKDIR) lib
	$(CXX) -shared -o lib/liblapack_wrapper_linux.so $(OBJS) $(LIBS)

lib/liblapack_wrapper_linux_static.a: $(OBJS)
	@$(MKDIR) lib
	$(AR) lib/liblapack_wrapper_linux_static.a $(OBJS)
