#ifndef __TYPEDEFINT
#define __TYPEDEFINT
 
#if defined (LONGBLAS) || defined (BLAS64) || defined (MATLAB_MEX_FILE)
typedef long Int;
#else
typedef int Int;
#endif
 
#endif