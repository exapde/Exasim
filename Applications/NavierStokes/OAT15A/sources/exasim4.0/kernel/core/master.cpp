#ifndef __MASTER
#define __MASTER

#include <vector>
#include <math.h>

using namespace std;
typedef int Int; 

extern "C" {
    void DGEMM(char*,char*,Int*,Int*,Int*,double*,double*,Int*,
             double*,Int*,double*,double*,Int*);        
    void DGETRF(Int*,Int*,double*,Int*,Int*,Int*);
    void DGETRI(Int*,double*,Int*,Int*,double*,Int*,Int*);
    void DSYEV(char*,char*,Int*,double*,Int*,double*,double*,Int*,Int*);
}

#include "quadrature.cpp"
#include "masternodes.cpp"
#include "mkshape.cpp"

#endif

