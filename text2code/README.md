# text2code

cd text2code 

mkdir build

cd build

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$HOME/GitHub/text2code \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=OFF \
  -DWITH_GMP=OFF \
  -DWITH_MPFR=OFF \
  -DWITH_LLVM=OFF \
  -DWITH_SYMENGINE_THREAD_SAFE=ON \
  -DINTEGER_CLASS=boostmp

make -j4

make install 

cd ../text2code 

g++ -O2 -std=c++17 Text2Code.cpp -o text2code

cd ../examples 

./../text2code/text2code pdemodel.txt


