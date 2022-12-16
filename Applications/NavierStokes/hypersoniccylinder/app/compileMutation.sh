#!/bin/bash
g++ -fPIC -O3 -c opuApp.cpp -w -std=c++11 -D _MUTATIONPP -I /Users/rloekvh/Mutationpp/install/include/mutation++/ -I /Users/rloekvh/Mutationpp/thirdparty/eigen/ -I /Users/rloekvh/Mutationpp/install
ar -rvs opuApp.a opuApp.o
g++ -std=c++11 ../../../../src/Kernel/Main/main.cpp -w -o serialapp ../../../../lib/Mac/commonCore.a ../../../../lib/Mac/opuCore.a opuApp.a -O2 -ldl -lm -lblas -llapack -lmutation++ -std=c++11 -D _MUTATIONPP -L /Users/rloekvh/Mutationpp/install/lib/ -I /Users/rloekvh/Mutationpp/install/include/mutation++/ -I /Users/rloekvh/Mutationpp/thirdparty/eigen/ -I /Users/rloekvh/Mutationpp/install/include/
