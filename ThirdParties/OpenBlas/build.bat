mkdir build2
cd build2
cmake -DUSE_THREAD:BOOL=false -DDYNAMIC_ARCH:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON ..
cmake --build . --clean-first --config Release --target install