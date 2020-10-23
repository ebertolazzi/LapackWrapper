BITS    = x64
WARN    = -Wall
CC      = gcc $(WARN)
CXX     = g++ -std=c++11 -pthread $(WARN)
AR      = ar rcs
LIBSGCC = -lstdc++ -lm -pthread

override LIBS += lib3rd/lib/libUtils_static_$(BITS).lib 

ifneq (,$(findstring LAPACK_WRAPPER_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring LAPACK_WRAPPER_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += lib3rd/lib/libopenblas_$(BITS).lib -Llib3rd/lib -Llib3rd/dll/$(BITS) -L$(FPATH) -lgfortran
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

ALL_LIBS = $(LIBS) -Llib/lib  -Llib/dll -llapack_wrapper_mingw_$(BITS).dll -lHSL_mingw_$(BITS).dll -lUtils_mingw$(BITS).dll
# -lHSL_mingw_$(BITS)

override DEFS += -DMINGW

all_libs: config lib/dll/liblapack_wrapper_mingw_$(BITS).dll.a lib/lib/liblapack_wrapper_mingw_$(BITS)_static.a lib/dll/libHSL_mingw_$(BITS).dll.a
#all_libs: config lib/liblapack_wrapper_mingw_$(BITS)_static.lib

lib/dll/liblapack_wrapper_mingw_$(BITS).dll.a: $(OBJS) lib/dll/libHSL_mingw_$(BITS).dll.a
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/bin
	@$(MKDIR) lib/dll
	$(CXX) -shared -o lib/bin/liblapack_wrapper_mingw_$(BITS).dll $(OBJS) \
	-Wl,--out-implib=lib/dll/liblapack_wrapper_mingw_$(BITS).dll.a \
	-Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--no-whole-archive \
	$(LIBS) lib/dll/libHSL_mingw_$(BITS).dll.a lib/dll/libUtils_mingw_$(BITS).dll.a

lib/lib/liblapack_wrapper_mingw_$(BITS)_static.a: $(OBJS)
	@$(MKDIR) lib
	$(AR) lib/lib/liblapack_wrapper_mingw_$(BITS)_static.a $(OBJS)

lib/dll/libHSL_mingw_$(BITS).dll.a: $(OBJS)
	@$(MKDIR) lib
	@$(MKDIR) lib/lib
	@$(MKDIR) lib/bin
	@$(MKDIR) lib/dll
	$(CXX) src/HSL/hsl_fake.cc -shared -o lib/bin/libHSL_mingw_$(BITS).dll \
	-Wl,--out-implib=lib/dll/libHSL_mingw_$(BITS).dll.a \
	-Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--no-whole-archive \
	$(LIBS)
