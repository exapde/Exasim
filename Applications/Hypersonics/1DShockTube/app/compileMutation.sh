#!/bin/bash
# /opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++
/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++ -fPIC -O3 -c opuApp.cpp -D _MUTATIONPP -I /Users/rloekvh/Mutationpp/install/include/mutation++/ -I /Users/rloekvh/Mutationpp/thirdparty/eigen/ -I /Users/rloekvh/Mutationpp/install/
ar -rvs opuApp.a opuApp.o
/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++ -std=c++11 ../../../../src/Kernel/Main/main.cpp -o serialapp ../../../../lib/Mac/commonCore.a ../../../../lib/Mac/opuCore.a opuApp.a -arch arm64 -O2 -ldl -lm -lblas -llapack -D _MUTATIONPP -I /Users/rloekvh/Mutationpp/install/include/mutation++/ -I /Users/rloekvh/Mutationpp/thirdparty/eigen/ -I /Users/rloekvh/Mutationpp/install/ -L /Users/rloekvh/Mutationpp/install/lib/ -lmutation++
