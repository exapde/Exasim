// backend/Common/common.cpp
//
// Definitions for the BLAS / CUBLAS / HIPBLAS sentinel globals declared
// `extern` in common.h. Carved out as part of the library port (see
// LIBRARY_PORT_INVENTORY.md, "ODR landmines") so that common.h can be
// included from multiple translation units in exasim_core.
//
// One TU in the runtime must compile this file; today that TU is the
// per-executable `main.cpp` source list, but it will move into
// `exasim_core` once that static library exists (Phase 1.3).
#include "common.h"

dstype one      =  1.0;
dstype minusone = -1.0;
dstype zero     =  0.0;

char chn = 'N';
char cht = 'T';
char chl = 'L';
char chu = 'U';
char chr = 'R';
char chv = 'V';

Int inc1 = 1;

dstype cublasOne     [1] = { 1.0};
dstype cublasMinusone[1] = {-1.0};
dstype cublasZero    [1] = { 0.0};
