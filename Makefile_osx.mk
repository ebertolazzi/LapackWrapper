WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
CC       = clang
CXX      = clang++
VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
#---------
CXX     += -std=c++11 -stdlib=libc++ -g
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
  override LIBS += -Llib3rd/lib -Llib3rd/dll -Wl,-rpath,lib3rd/dll -Wl,-rpath,lib/lib -L$(FPATH)/../../..  -lgfortran
  override LIBS += -L/usr/local/opt/openblas/lib -lopenblas
  override INC  += -I/usr/local/opt/openblas/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ATLAS,$(USED_LIB)))
  $(error error is "Atlas is supported only on Linux!")
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_MKL,$(USED_LIB)))
  # for MKL
  INTEL_PATH = /opt/intel/lib
  MKL_PATH   = /opt/intel/mkl
  MKL_LIB = -lmkl_sequential -lmkl_intel_thread -lmkl_core -liomp5
  override LIBS += -L$(MKL_PATH)/lib -Xlinker -rpath -Xlinker $(MKL_PATH)/lib -L$(INTEL_PATH) -Xlinker -rpath -Xlinker $(INTEL_PATH) $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_BLASFEO,$(USED_LIB)))
  override LIBS += -L/opt/blasfeo/lib -Wl,-rpath,/opt/blasfeo/lib -lblasfeo
  override INC  += -I/opt/blasfeo/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ACCELERATE,$(USED_LIB)))
  override LIBS += -framework Accelerate
endif

ALL_LIBS = $(LIBS) -Llib/lib -Llib/dll -llapack_wrapper_osx_static -lHSL_osx

all_libs: config lib/dll/liblapack_wrapper_osx.dylib lib/lib/liblapack_wrapper_osx_static.a lib/dll/libHSL_osx.dylib
##all_libs: config lib/liblapack_wrapper_osx_static.a

lib/dll/liblapack_wrapper_osx.dylib: lib/dll/libHSL_osx.dylib $(OBJS)
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	@$(MKDIR) lib/bin
	$(CXX) -dynamiclib -o lib/dll/liblapack_wrapper_osx.dylib lib/dll/libHSL_osx.dylib $(OBJS) $(LIBS) -Llib/lib -Llib/dll -lHSL_osx

lib/lib/liblapack_wrapper_osx_static.a: $(OBJS)
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	@$(MKDIR) lib/bin
	libtool -static $(OBJS) -o lib/lib/liblapack_wrapper_osx_static.a

lib/dll/libHSL_osx.dylib:
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	@$(MKDIR) lib/bin
	$(CXX) -dynamiclib -o lib/dll/libHSL_osx.dylib src/HSL/hsl_fake.cc $(LIBS)
