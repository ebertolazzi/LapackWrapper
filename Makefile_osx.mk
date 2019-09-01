WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
CC       = clang
CXX      = clang++
VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
#---------
CXX     += -std=c++11 -stdlib=libc++
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
  override LIBS += -L/usr/local/opt/openblas/lib
  override INC  += -I/usr/local/opt/openblas/include
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
  override LIBS += -framework Accelerate
endif

ALL_LIBS = $(LIBS) -Llib -llapack_wrapper_osx_static

##all_libs: lib/liblapack_wrapper_osx.dylib lib/liblapack_wrapper_osx_static.a
all_libs: config lib/liblapack_wrapper_osx_static.a

lib/liblapack_wrapper_osx.dylib: $(OBJS)
	@$(MKDIR) lib
	$(CXX) -shared -o lib/liblapack_wrapper_osx.dylib $(OBJS) $(LIBS)

lib/liblapack_wrapper_osx_static.a: $(OBJS)
	@$(MKDIR) lib
	libtool -static $(OBJS) -o lib/liblapack_wrapper_osx_static.a
