#ifndef __PREPROCESSING
#define __PREPROCESSING

#include "../preprocessing/errormsg.cpp"
#include "../preprocessing/ioutilities.cpp"
#include "../preprocessing/readbinaryfiles.cpp"

#ifdef HAVE_CUDA
#include "../preprocessing/gpuDeviceInfo.cpp"
#endif

#endif