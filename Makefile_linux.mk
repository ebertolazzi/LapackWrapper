WARN   = -Wall
CC     = gcc -fPIC
CXX    = g++ -fPIC -std=c++11 -pthread -g
#
# activate C++11 for g++ >= 4.9
#
VERSION  = $(shell $(CC) -dumpversion)
CC     += $(WARN)
CXX    += $(WARN)
AR      = ar rcs
LIBSGCC = -lstdc++ -lm -pthread

override LIBS += -lfort_linux_static -lfmt_linux_static

ifneq (,$(findstring LAPACK_WRAPPER_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += -Llib3rd/dll -Llib3rd/lib -Wl,-rpath,./lib/dll -lopenblas -L$(FPATH)/../../.. -lgfortran
  override INC  += -I/usr/include/x86_64-linux-gnu/
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

ifneq (,$(findstring LAPACK_WRAPPER_USE_BLASFEO,$(USED_LIB)))
  override LIBS += -L/opt/blasfeo/lib -Wl,-rpath,/opt/blasfeo/lib -lblasfeo
  override INC  += -I/opt/blasfeo/include
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_ACCELERATE,$(USED_LIB)))
  $(error error is "Accelerate is supported only on Darwin!")
endif

ALL_LIBS = $(LIBS) -Llib/lib -Llib/dll -llapack_wrapper_linux -lHSL_linux

all_libs: config lib/dll/liblapack_wrapper_linux.so lib/lib/liblapack_wrapper_linux_static.a lib/dll/libHSL_linux.so

lib/dll/liblapack_wrapper_linux.so: lib/dll/libHSL_linux.so $(OBJS)
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	$(CXX) -shared -o lib/dll/liblapack_wrapper_linux.so $(OBJS) $(LIBS) -Llib/dll -Llib/lib -lHSL_linux

lib/lib/liblapack_wrapper_linux_static.a: $(OBJS)
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	$(AR) lib/lib/liblapack_wrapper_linux_static.a $(OBJS)

lib/dll/libHSL_linux.so:
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/dll
	$(CXX) -shared -o lib/dll/libHSL_linux.so src/HSL/hsl_fake.cc $(LIBS)
