mkdir build2
cd build2
cmake -DUSE_THREAD:VAR=false -DDYNAMIC_ARCH:VAR=ON -DBUILD_SHARED_LIBS:VAR=ON ..
cmake --build . --clean-first --config Release --target install