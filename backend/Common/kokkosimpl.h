#ifndef __KOKKOSIMPL_H__
#define __KOKKOSIMPL_H__

void AverageFlux(dstype* fg, const int N)
{	
    Kokkos::parallel_for("AverageFlux", N, KOKKOS_LAMBDA(const size_t tid) {
        fg[tid+N] = 0.5*(fg[tid] + fg[tid+N]);            
    });
}

void AverageFluxDotNormal(dstype* fg, const dstype* nl, const int N, const int M, const int numPoints, const int nd)
{	
    Kokkos::parallel_for("AverageFluxDotNormal", M, KOKKOS_LAMBDA(const size_t tid) {
        int i = tid%numPoints;                
        dstype sum = fg[N + 0*M + tid] * nl[i + 0 * numPoints];   
        for (int l = 1; l < nd; l++)
            sum += fg[N + l*M + tid] * nl[i + l * numPoints];
        fg[tid] = sum;     
    });
}

void FluxDotNormal(dstype* fh, dstype* fg, const dstype* nl, const int M, const int numPoints, const int nd)
{	
    Kokkos::parallel_for("FluxDotNormal", M, KOKKOS_LAMBDA(const size_t tid) {
        int i = tid%numPoints;                
        dstype sum = fg[0*M + tid] * nl[i + 0 * numPoints];   
        for (int l = 1; l < nd; l++)
            sum += fg[l*M + tid] * nl[i + l * numPoints];
        fh[tid] = sum;     
    });
}

// void FluxDerivativeDotNormal(dstype* fh, dstype* fg, const dstype* nl, const int M, const int numPoints, const int nd, const int nc)
// {	
//     Kokkos::parallel_for("FluxDotNormal", M, KOKKOS_LAMBDA(const size_t tid) {
//         int i = tid%numPoints;
//         for (int n = 0; n < nc; n++) {
//           dstype sum = fg[tid + M*0 + M*nd*n] * nl[i + 0 * numPoints];   
//           for (int l = 1; l < nd; l++)
//               sum += fg[tid + M*l + M*nd*n] * nl[i + l * numPoints];
//           fh[tid + M*n] = sum;     
//         }
//     });
// }

void AddStabilization1(dstype* fg, const dstype* ug1, const dstype* ug2, const dstype* tau, const int M)
{	
    Kokkos::parallel_for("AddStabilization1", M, KOKKOS_LAMBDA(const size_t tid) {
        fg[tid] += tau[0] * (ug1[tid] - ug2[tid]);        
    });
}

void AddStabilization1(dstype* fg, dstype* fg1, dstype* fg2, const dstype* ug1, const dstype* ug2, const dstype* tau, const int M, const int numPoints)
{	
    Kokkos::parallel_for("AddStabilization1", M, KOKKOS_LAMBDA(const size_t tid) {
        int i = tid%numPoints;
        int n = tid/numPoints;
        fg[tid] += tau[0] * (ug1[tid] - ug2[tid]);        
        fg1[i + numPoints*n + M*n] += tau[0]; 
        fg2[i + numPoints*n + M*n] -= tau[0]; 
    });
}

void AddStabilization2(dstype* fg, const dstype* ug1, const dstype* ug2, const dstype* tau, const int M, const int numPoints)
{	
    Kokkos::parallel_for("AddStabilization2", M, KOKKOS_LAMBDA(const size_t tid) {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
        fg[tid] += tau[j] * (ug1[tid] - ug2[tid]);     
    });
}

void AddStabilization3(dstype* fg, const dstype* ug1, const dstype* ug2, const dstype* tau, const int M, const int numPoints, const int ncu)
{	
    Kokkos::parallel_for("AddStabilization3", M, KOKKOS_LAMBDA(const size_t tid) {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
        for (int k=0; k<ncu; k++) {
            int nm = k * ncu + j;
            int nk = k * numPoints + i;
            fg[tid] += tau[nm] * (ug1[nk] - ug2[nk]);
        }
    });
}

void GetArrayAtIndex(dstype* y, const dstype* x, const int* ind, const int n)
{    
    Kokkos::parallel_for("GetArrayAtIndex", n, KOKKOS_LAMBDA(const size_t i) {
        y[i] = x[ind[i]];
    });
}

void PutArrayAtIndex(dstype* y, const dstype* x, const int* ind, const int n)
{    
    Kokkos::parallel_for("PutArrayAtIndex", n, KOKKOS_LAMBDA(const size_t i) {
        y[ind[i]] = x[i];
    });
}

void AddColumns(dstype* y, const dstype* x, const int m, const int n)
{    
    Kokkos::parallel_for("GetCollumnAtIndex", m, KOKKOS_LAMBDA(const size_t idx) {        
      for (int j=0; j<n; j++) y[idx] += x[idx + m*j];
    });
}

void SubtractColumns(dstype* y, const dstype* x, const int m, const int n)
{    
    Kokkos::parallel_for("GetCollumnAtIndex", m, KOKKOS_LAMBDA(const size_t idx) {        
      for (int j=0; j<n; j++) y[idx] -= x[idx + m*j];
    });
}

void GetCollumnAtIndex(dstype* y, const dstype* x, const int* ind, const int m, const int n)
{    
    int N = m*n;
    Kokkos::parallel_for("GetCollumnAtIndex", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%m;
        int j = idx/m;
        y[idx] = x[i + m*ind[j]];
    });
}

void PutCollumnAtIndex(dstype* y, const dstype* x, const int* ind, const int m, const int n)
{    
    int N = m*n;
    Kokkos::parallel_for("GetCollumnAtIndex", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%m;
        int j = idx/m;
        y[i + m*ind[j]] = x[idx];
    });
}

void PutCollumnAtIndexAtomicAdd(dstype* y, const dstype* x, const int* ind, const int m, const int n)
{    
    int N = m*n;
    Kokkos::parallel_for("GetCollumnAtIndex", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%m;
        int j = idx/m;
        //y[i + m*ind[j]] = x[idx];
        Kokkos::atomic_add(&y[i + m*ind[j]], x[idx]);   
    });
}

void GetCollumnAtIndex(dstype* y, const dstype* x, const int* ind, const int i0, const int k, const int m, const int n)
{    
    int N = m*n;
    Kokkos::parallel_for("GetCollumnAtIndex", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%m;
        int j = idx/m;
        y[i0 + i + k*j] = x[i + m*ind[j]];
    });
}

void PutCollumnAtIndex(dstype* y, const dstype* x, const int* ind, const int i0, const int k, const int m, const int n)
{    
    int N = m*n;
    Kokkos::parallel_for("GetCollumnAtIndex", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%m;
        int j = idx/m;
        y[i + m*ind[j]] = x[i0 + i + k*j];
    });
}

void ArrayCopy(dstype* y, const dstype* x, const int n)
{    
    Kokkos::parallel_for("ArrayCopy", n, KOKKOS_LAMBDA(const size_t i) {
        y[i] = x[i];
    });
}

void ArraySetValue(dstype* y, const dstype a, const int n)
{    
    Kokkos::parallel_for("ArraySetValue", n, KOKKOS_LAMBDA(const size_t i) {
        y[i] = a;
    });
}

void ArrayAddScalar(dstype* y, const dstype a, const int n)
{    
    Kokkos::parallel_for("ArraySetValue", n, KOKKOS_LAMBDA(const size_t i) {
        y[i] += a;
    });
}

void ArrayMultiplyScalar(dstype* y, const dstype a, const int n)
{    
    Kokkos::parallel_for("ArrayMultiplyScalar", n, KOKKOS_LAMBDA(const size_t i) {
        y[i] = a*y[i];
    });
}

void ArrayAXPB(dstype* y, dstype* x, const dstype a, const dstype b, const int N) 
{
    Kokkos::parallel_for("ArrayAXPBY", N, KOKKOS_LAMBDA(const size_t idx) {
        y[idx] = a * x[idx] + b;
    });
}

void ArrayAXPBY(dstype* y, dstype* x, const dstype* z, const dstype a, const dstype b, const int N) 
{
    Kokkos::parallel_for("ArrayAXPBY", N, KOKKOS_LAMBDA(const size_t idx) {
        y[idx] = a * x[idx] + b * z[idx];
    });
}

void ArrayAXY(dstype* y, dstype* x, const dstype* z, const dstype a, const int N) 
{
    Kokkos::parallel_for("ArrayAXY", N, KOKKOS_LAMBDA(const size_t idx) {
        y[idx] = a * x[idx] * z[idx];
    });
}

void ArrayAdd3Vectors(dstype* s, dstype* x, dstype* y, dstype* z, dstype a, dstype b, dstype c, int N)
{    
    Kokkos::parallel_for("ArrayAXY", N, KOKKOS_LAMBDA(const size_t idx) {
        s[idx] = a*x[idx] + b*y[idx] + c*z[idx];   
    });
}

void ArrayExtract(dstype* un, const dstype* u, const int I, const int J, const int K, 
        const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;
    Kokkos::parallel_for("ArrayExtract", N, KOKKOS_LAMBDA(const size_t idx) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+Q*k];
    });
}

void ArrayInsert(dstype* u, const dstype* un, const int I, const int J, const int K, 
        const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;
    Kokkos::parallel_for("ArrayInsert", N, KOKKOS_LAMBDA(const size_t idx) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+Q*k] = un[idx];        
    });
}


void UpdateUDG(dstype* u, const dstype* un, const dstype alpha, const int I, const int J, const int K, 
        const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;
    Kokkos::parallel_for("UpdateUDG", N, KOKKOS_LAMBDA(const size_t idx) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+Q*k] += alpha*un[idx];        
    });
}

void columnwiseMultiply(dstype* C, const dstype* A, const dstype* b, const int N, const int M)
{    
    int K = N*M;
    Kokkos::parallel_for("columnwiseMultiply", K, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%N;   // [1, N]         
        C[idx] = A[idx]*b[i];           
    });    
}

void ArrayGemmBatch(dstype* C, const dstype* A, const dstype* B, const int I, const int J, const int K, const int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    Kokkos::parallel_for("ArrayGemmBatch", N, KOKKOS_LAMBDA(const size_t idx) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        dstype sum = 0.0;
        for (int k=0; k<K; k++)
            sum += A[i+I*k+I*K*s]*B[k+K*j+K*J*s];
        C[idx] = sum;    
    });
}

void ArrayGemmBatch1(dstype* C, const dstype* A, const dstype* B, const int I, const int J, const int K, const int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    Kokkos::parallel_for("ArrayGemmBatch1", N, KOKKOS_LAMBDA(const size_t idx) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        dstype sum = C[idx];
        for (int k=0; k<K; k++)
            sum += A[i+I*k+Q*s]*B[k+K*j+P*s];
       C[idx] = sum;     
    });
}

void ArrayGemmBatch2(dstype* C, const dstype* A, const dstype* B, dstype alpha, const int I, const int J, const int K, const int S)
{        
    // C[S*I*J] = A[S*I*K] x B[S*K*J] + C[S*I*J]
    int M = I*J;
    int N = M*S;
    int Q = S*I;
    int P = S*K;
    Kokkos::parallel_for("ArrayGemmBatch2", N, KOKKOS_LAMBDA(const size_t idx) {
        int s = idx%S;        
        int l = idx/S;
        int i = l%I;
        int j = l/I;
        dstype sum = alpha*C[idx];
        for (int k=0; k<K; k++)
            sum += A[s+S*i+Q*k]*B[s+S*k+P*j];
        C[idx] = sum;     
    });
}

void ArrayDG2CG(dstype* ucg, const dstype* udg, const int* cgent2dgent, const int* rowent2elem, const int nent)
{        
    Kokkos::parallel_for("ArrayDG2CG", nent, KOKKOS_LAMBDA(const size_t i) {
        dstype sum = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        dstype fac = 1.0/((dstype) nelem);
        for (int k=0; k<nelem; k++)
            sum += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = sum*fac;
    });
}

void ArrayDG2CG2(dstype* ucg, const dstype* udg, const int* colent2elem, const int* rowent2elem, const int nent, const int npe)
{        
    Kokkos::parallel_for("ArrayDG2CG2", nent, KOKKOS_LAMBDA(const size_t i) {        
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        dstype fac = 1.0/((dstype) (nelem*npe));
        dstype sum = 0.0;
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            for (int j=0; j<npe; j++)
                sum += udg[j+npe*e]; 
        }
        ucg[i] = sum*fac;
    });
}

void computeQTv(dstype* p, const dstype* Q, const dstype* v, const int M, const int N)
{
    // Compute p = Q^T v
    Kokkos::parallel_for("ComputeQTv", Kokkos::TeamPolicy<>(N, Kokkos::AUTO),
        KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
        int col = team.league_rank();
        double sum = 0.0;

        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team, M),
            [&](const int i, double& local_sum) {
                local_sum += Q[i + M * col] * v[i];
            }, sum);

        if (team.team_rank() == 0) {
            p[col] = sum;
        }
    });
}

void ArrayMatrixMultiplication(dstype* C, const dstype* A, const dstype* B, const int I, const int J, const int K)
{        
    // C[I*J] = A[I*K] x B[K*J]
    int N = I*J;    
    Kokkos::parallel_for("ArrayMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%I;   //   [1, I]
        int j = idx/I; // [1, J]        
        dstype sum = 0.0;
        for (int k=0; k<K; k++)
            sum += A[i + I*k]*B[k + K*j];
        C[idx] = sum;    
    });
}

void ArrayMatrixMultiplication1(dstype* C, const dstype* A, const dstype* B, const int I, const int J, const int K)
{        
    // C[I*J] += A[I*K] x B[K*J]
    int N = I*J;    
    Kokkos::parallel_for("ArrayMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%I;   //   [1, I]
        int j = idx/I; // [1, J]        
        dstype sum = 0.0;
        for (int k=0; k<K; k++)
            sum += A[i + I*k]*B[k + K*j];
        C[idx] += sum;    
    });
}

void ArrayEosInverseMatrix11(dstype* A, const int npe, const int ncw, const int ne)
{            
    int N = npe*ne;
    Kokkos::parallel_for("ArrayEosInverseMatrix11", N, KOKKOS_LAMBDA(const size_t i) {
        A[i] = 1.0/A[i];    
    });
}

// void AVdistfunc(dstype* A, const dstype* param, const int npe, const int ncu, const int ne)
// {        
//     int N = npe*ne;
//     Kokkos::parallel_for("AVdistfunc", N, KOKKOS_LAMBDA(const size_t i) {
//         int j = i%npe; // [1, npe]
//         int k = i/npe; // [1, ne]
//         A[j + npe*ncu*k] = param[0] * tanh(param[1] * A[j + npe + npe*ncu*k]);
//     });
// }

void ArrayEosInverseMatrix22(dstype* A, const int npe, const int ncw, const int ne)
{        
    int N = npe*ne;
    int M = npe*ncw*ncw;
    Kokkos::parallel_for("ArrayEosInverseMatrix22", N, KOKKOS_LAMBDA(const size_t i) {
        int j = i%npe;    // [1, npe]
        int k = (i-j)/npe; //[1, ne]
        int jk = j + M*k;
        dstype a11 = A[jk + npe*0];
        dstype a21 = A[jk + npe*1];
        dstype a12 = A[jk + npe*2];
        dstype a22 = A[jk + npe*3];
        dstype detA = (a11*a22- a12*a21);      
        A[jk + npe*0] = a22/detA;
        A[jk + npe*1] = -a21/detA;
        A[jk + npe*2] = -a12/detA;
        A[jk + npe*3] = a11/detA;
    });
}

void ArrayEosInverseMatrix33(dstype* A, const int npe, const int ncw, const int ne)
{            
    int N = npe*ne;
    int M = npe*ncw*ncw;
    Kokkos::parallel_for("ArrayEosInverseMatrix33", N, KOKKOS_LAMBDA(const size_t i) {
        int j = i%npe;    // [1, npe]
        int k = (i-j)/npe; //[1, ne]
        int jk = j + M*k;
        dstype a11 = A[jk + npe*0];
        dstype a21 = A[jk + npe*1];
        dstype a31 = A[jk + npe*2];
        dstype a12 = A[jk + npe*3];
        dstype a22 = A[jk + npe*4];
        dstype a32 = A[jk + npe*5];
        dstype a13 = A[jk + npe*6];
        dstype a23 = A[jk + npe*7];
        dstype a33 = A[jk + npe*8];        
        dstype detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
      
        A[jk + npe*0] = (a22*a33 - a23*a32)/detA;
        A[jk + npe*1] = (a23*a31 - a21*a33)/detA;
        A[jk + npe*2] = (a21*a32 - a22*a31)/detA;
        A[jk + npe*3] = (a13*a32 - a12*a33)/detA;
        A[jk + npe*4] = (a11*a33 - a13*a31)/detA;
        A[jk + npe*5] = (a12*a31 - a11*a32)/detA;
        A[jk + npe*6] = (a12*a23 - a13*a22)/detA;
        A[jk + npe*7] = (a13*a21 - a11*a23)/detA;
        A[jk + npe*8] = (a11*a22 - a12*a21)/detA;
    });
}

void ArrayEosMatrixMultiplication(dstype* C, const dstype* A, const dstype* B, const int npe, const int ncw, const int ne, const int ncu)
{        
    // C[npe*ncw*ncu*ne] = A[npe*ncw*ncw*ne] x B[npe*ncw*ncu*ne]
    int N = npe*ne;
    int K = npe*ncw;
    int P = K*ncu;
    int Q = K*ncw;    
    Kokkos::parallel_for("ArrayEosMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t i) {
        int j = i%npe;    // [1, npe]
        int k = (i-j)/npe; //[1, ne]        
        for (int b=0; b<ncu; b++)
          for (int a=0; a<ncw; a++) {
            dstype sum = 0.0;
            for (int m=0; m<ncw; m++)
              sum += A[j + npe*a + K*m + Q*k]*B[j + npe*m + K*b + P*k];
            C[j + npe*a + K*b + P*k] = sum;  
          }        
    });
}

void SmallMatrixSolve11(dstype *b, dstype* A, const int N, const int nrhs)
{                
    Kokkos::parallel_for("SmallMatrixSolve11", N, KOKKOS_LAMBDA(const size_t i) {        
        dstype c1 = 1.0/A[i];
        for (int j=0; j<nrhs; j++)
            b[i + N*j] = c1*b[i + N*j];        
    });
}


void SmallMatrixSolve22(dstype *b, dstype* A, const int N, const int nrhs)
{                
    Kokkos::parallel_for("SmallMatrixSolve22", N, KOKKOS_LAMBDA(const size_t i) {        
        dstype a11 = A[i + N*0];
        dstype a21 = A[i + N*1];
        dstype a12 = A[i + N*2];
        dstype a22 = A[i + N*3];
        dstype detA = (a11*a22- a12*a21);              
        dstype c11 = a22/detA;
        dstype c21 = -a21/detA;
        dstype c12 = -a12/detA;
        dstype c22 = a11/detA;
        for (int j=0; j<nrhs; j++) {
          dstype b1 = b[i + N*0 + N*2*j];
          dstype b2 = b[i + N*1 + N*2*j];        
          b[i + N*0 + N*2*j] = c11*b1 + c12*b2;
          b[i + N*1 + N*2*j] = c21*b1 + c22*b2;
        }        
    });
}

void SmallMatrixSolve33(dstype *b, dstype* A, const int N, const int nrhs)
{                
    Kokkos::parallel_for("SmallMatrixSolve22", N, KOKKOS_LAMBDA(const size_t i) {      

        dstype a11 = A[i + N*0];
        dstype a21 = A[i + N*1];
        dstype a31 = A[i + N*2];
        dstype a12 = A[i + N*3];
        dstype a22 = A[i + N*4];
        dstype a32 = A[i + N*5];
        dstype a13 = A[i + N*6];
        dstype a23 = A[i + N*7];
        dstype a33 = A[i + N*8];        
        dstype detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);

        dstype c11 = (a22*a33 - a23*a32)/detA;
        dstype c21 = (a23*a31 - a21*a33)/detA;
        dstype c31 = (a21*a32 - a22*a31)/detA;
        dstype c12 = (a13*a32 - a12*a33)/detA;
        dstype c22 = (a11*a33 - a13*a31)/detA;
        dstype c32 = (a12*a31 - a11*a32)/detA;
        dstype c13 = (a12*a23 - a13*a22)/detA;
        dstype c23 = (a13*a21 - a11*a23)/detA;
        dstype c33 = (a11*a22 - a12*a21)/detA;

        if (nrhs == 1) {
          dstype b1 = b[i + N*0];
          dstype b2 = b[i + N*1];
          dstype b3 = b[i + N*2];

          b[i + N*0] = c11*b1 + c12*b2 + c13*b3;
          b[i + N*1] = c21*b1 + c22*b2 + c23*b3;
          b[i + N*2] = c31*b1 + c32*b2 + c33*b3;        
        }
        else {
          for (int j=0; j<nrhs; j++) {
            dstype b1 = b[i + N*0 + N*3*j];
            dstype b2 = b[i + N*1 + N*3*j];     
            dstype b3 = b[i + N*2 + N*3*j];   
            b[i + N*0 + N*3*j] = c11*b1 + c12*b2 + c13*b3;
            b[i + N*1 + N*3*j] = c21*b1 + c22*b2 + c23*b3;
            b[i + N*2 + N*3*j] = c31*b1 + c32*b2 + c33*b3;
          }        
        }
    });
}

void SmallMatrixSolve44(dstype *b, dstype* A, const int N, const int nrhs)
{                
  Kokkos::parallel_for("SmallMatrixSolve44", N, KOKKOS_LAMBDA(const size_t i) {      

    dstype a11 = A[i + N*0];
    dstype a21 = A[i + N*1];
    dstype a31 = A[i + N*2];
    dstype a41 = A[i + N*3];
    dstype a12 = A[i + N*4];
    dstype a22 = A[i + N*5];
    dstype a32 = A[i + N*6];
    dstype a42 = A[i + N*7];
    dstype a13 = A[i + N*8];
    dstype a23 = A[i + N*9];
    dstype a33 = A[i + N*10];
    dstype a43 = A[i + N*11];
    dstype a14 = A[i + N*12];
    dstype a24 = A[i + N*13];
    dstype a34 = A[i + N*14];
    dstype a44 = A[i + N*15];        
    dstype detA = (a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42
             - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41
             + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41
             - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

    dstype c11 = (a22*a33*a44 - a22*a34*a43 - a23*a32*a44 + a23*a34*a42 + a24*a32*a43 - a24*a33*a42) / detA;
    dstype c21 = (a21*a33*a44 - a21*a34*a43 - a23*a31*a44 + a23*a34*a41 + a24*a31*a43 - a24*a33*a41) / detA;
    dstype c31 = (a21*a32*a44 - a21*a34*a42 - a22*a31*a44 + a22*a34*a41 + a24*a31*a42 - a24*a32*a41) / detA;
    dstype c41 = (a21*a32*a43 - a21*a33*a42 - a22*a31*a43 + a22*a33*a41 + a23*a31*a42 - a23*a32*a41) / detA;
    dstype c12 = (a12*a33*a44 - a12*a34*a43 - a13*a32*a44 + a13*a34*a42 + a14*a32*a43 - a14*a33*a42) / detA;
    dstype c22 = (a11*a33*a44 - a11*a34*a43 - a13*a31*a44 + a13*a34*a41 + a14*a31*a43 - a14*a33*a41) / detA;
    dstype c32 = (a11*a32*a44 - a11*a34*a42 - a12*a31*a44 + a12*a34*a41 + a14*a31*a42 - a14*a32*a41) / detA;
    dstype c42 = (a11*a32*a43 - a11*a33*a42 - a12*a31*a43 + a12*a33*a41 + a13*a31*a42 - a13*a32*a41) / detA;
    dstype c13 = (a12*a23*a44 - a12*a24*a43 - a13*a22*a44 + a13*a24*a42 + a14*a22*a43 - a14*a23*a42) / detA;
    dstype c23 = (a11*a23*a44 - a11*a24*a43 - a13*a21*a44 + a13*a24*a41 + a14*a21*a43 - a14*a23*a41) / detA;
    dstype c33 = (a11*a22*a44 - a11*a24*a42 - a12*a21*a44 + a12*a24*a41 + a14*a21*a42 - a14*a22*a41) / detA;
    dstype c43 = (a11*a22*a43 - a11*a23*a42 - a12*a21*a43 + a12*a23*a41 + a13*a21*a42 - a13*a22*a41) / detA;
    dstype c14 = (a12*a23*a34 - a12*a24*a33 - a13*a22*a34 + a13*a24*a32 + a14*a22*a33 - a14*a23*a32) / detA;
    dstype c24 = (a11*a23*a34 - a11*a24*a33 - a13*a21*a34 + a13*a24*a31 + a14*a21*a33 - a14*a23*a31) / detA;
    dstype c34 = (a11*a22*a34 - a11*a24*a32 - a12*a21*a34 + a12*a24*a31 + a14*a21*a32 - a14*a22*a31) / detA;
    dstype c44 = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31) / detA;

    if (nrhs == 1) {
      dstype b1 = b[i + N*0];
      dstype b2 = b[i + N*1];
      dstype b3 = b[i + N*2];
      dstype b4 = b[i + N*3];

      b[i + N*0] = c11*b1 + c12*b2 + c13*b3 + c14*b4;
      b[i + N*1] = c21*b1 + c22*b2 + c23*b3 + c24*b4;
      b[i + N*2] = c31*b1 + c32*b2 + c33*b3 + c34*b4;
      b[i + N*3] = c41*b1 + c42*b2 + c43*b3 + c44*b4;        
    }
    else {
      for (int j=0; j<nrhs; j++) {
        dstype b1 = b[i + N*0 + N*4*j];
        dstype b2 = b[i + N*1 + N*4*j];     
        dstype b3 = b[i + N*2 + N*4*j];   
        dstype b4 = b[i + N*3 + N*4*j];
        b[i + N*0 + N*4*j] = c11*b1 + c12*b2 + c13*b3 + c14*b4;
        b[i + N*1 + N*4*j] = c21*b1 + c22*b2 + c23*b3 + c24*b4;
        b[i + N*2 + N*4*j] = c31*b1 + c32*b2 + c33*b3 + c34*b4;
        b[i + N*3 + N*4*j] = c41*b1 + c42*b2 + c43*b3 + c44*b4;
      }              
    }
  });
}

void SmallMatrixSolve55(dstype *b, dstype* A, const int N, const int nrhs)
{                
  Kokkos::parallel_for("SmallMatrixSolve55", N, KOKKOS_LAMBDA(const size_t i) {      

    dstype a11 = A[i + N*0];
    dstype a21 = A[i + N*1];
    dstype a31 = A[i + N*2];
    dstype a41 = A[i + N*3];
    dstype a51 = A[i + N*4];
    dstype a12 = A[i + N*5];
    dstype a22 = A[i + N*6];
    dstype a32 = A[i + N*7];
    dstype a42 = A[i + N*8];
    dstype a52 = A[i + N*9];
    dstype a13 = A[i + N*10];
    dstype a23 = A[i + N*11];
    dstype a33 = A[i + N*12];
    dstype a43 = A[i + N*13];
    dstype a53 = A[i + N*14];
    dstype a14 = A[i + N*15];
    dstype a24 = A[i + N*16];
    dstype a34 = A[i + N*17];
    dstype a44 = A[i + N*18];
    dstype a54 = A[i + N*19];
    dstype a15 = A[i + N*20];
    dstype a25 = A[i + N*21];
    dstype a35 = A[i + N*22];
    dstype a45 = A[i + N*23];
    dstype a55 = A[i + N*24];        
    dstype detA = (a11*a22*a33*a44*a55 - a11*a22*a33*a45*a54 - a11*a22*a34*a43*a55 + a11*a22*a34*a45*a53 + a11*a22*a35*a43*a54 - a11*a22*a35*a44*a53
             - a11*a23*a32*a44*a55 + a11*a23*a32*a45*a54 + a11*a23*a34*a42*a55 - a11*a23*a34*a45*a52 - a11*a23*a35*a42*a54 + a11*a23*a35*a44*a52
             + a11*a24*a32*a43*a55 - a11*a24*a32*a45*a53 - a11*a24*a33*a42*a55 + a11*a24*a33*a45*a52 + a11*a24*a35*a42*a53 - a11*a24*a35*a43*a52
             - a11*a25*a32*a43*a54 + a11*a25*a32*a44*a53 + a11*a25*a33*a42*a54 - a11*a25*a33*a44*a52 - a11*a25*a34*a42*a53 + a11*a25*a34*a43*a52
             - a12*a21*a33*a44*a55 + a12*a21*a33*a45*a54 + a12*a21*a34*a43*a55 - a12*a21*a34*a45*a53 - a12*a21*a35*a43*a54 + a12*a21*a35*a44*a53
             + a12*a23*a31*a44*a55 - a12*a23*a31*a45*a54 - a12*a23*a34*a41*a55 + a12*a23*a34*a45*a51 + a12*a23*a35*a41*a54 - a12*a23*a35*a44*a51
             - a12*a24*a31*a43*a55 + a12*a24*a31*a45*a53 + a12*a24*a33*a41*a55 - a12*a24*a33*a45*a51 - a12*a24*a35*a41*a53 + a12*a24*a35*a43*a51
             + a12*a25*a31*a43*a54 - a12*a25*a31*a44*a53 - a12*a25*a33*a41*a54 + a12*a25*a33*a44*a51 + a12*a25*a34*a41*a53 - a12*a25*a34*a43*a51
             + a13*a21*a32*a44*a55 - a13*a21*a32*a45*a54 - a13*a21*a34*a42*a55 + a13*a21*a34*a45*a52 + a13*a21*a35*a42*a54 - a13*a21*a35*a44*a52
             - a13*a22*a31*a44*a55 + a13*a22*a31*a45*a54 + a13*a22*a34*a41*a55 - a13*a22*a34*a45*a51 - a13*a22*a35*a41*a54 + a13*a22*a35*a44*a51
             + a13*a24*a31*a42*a55 - a13*a24*a31*a45*a52 - a13*a24*a32*a41*a55 + a13*a24*a32*a45*a51 + a13*a24*a35*a41*a52 - a13*a24*a35*a42*a51
             - a13*a25*a31*a42*a54 + a13*a25*a31*a44*a52 + a13*a25*a32*a41*a54 - a13*a25*a32*a44*a51 - a13*a25*a34*a41*a52 + a13*a25*a34*a42*a51
             - a14*a21*a32*a43*a55 + a14*a21*a32*a45*a53 + a14*a21*a33*a42*a55 - a14*a21*a33*a45*a52 - a14*a21*a35*a42*a53 + a14*a21*a35*a43*a52
             + a14*a22*a31*a43*a55 - a14*a22*a31*a45*a53 - a14*a22*a33*a41*a55 + a14*a22*a33*a45*a51 + a14*a22*a35*a41*a53 - a14*a22*a35*a43*a51
             - a14*a23*a31*a42*a55 + a14*a23*a31*a45*a52 + a14*a23*a32*a41*a55 - a14*a23*a32*a45*a51 - a14*a23*a35*a41*a52 + a14*a23*a35*a42*a51
             + a14*a25*a31*a42*a53 - a14*a25*a31*a43*a52 - a14*a25*a32*a41*a53 + a14*a25*a32*a43*a51 + a14*a25*a33*a41*a52 - a14*a25*a33*a42*a51
             + a15*a21*a32*a43*a54 - a15*a21*a32*a44*a53 - a15*a21*a33*a42*a54 + a15*a21*a33*a44*a52 + a15*a21*a34*a42*a53 - a15*a21*a34*a43*a52
             - a15*a22*a31*a43*a54 + a15*a22*a31*a44*a53 + a15*a22*a33*a41*a54 - a15*a22*a33*a44*a51 - a15*a22*a34*a41*a53 + a15*a22*a34*a43*a51
             + a15*a23*a31*a42*a54 - a15*a23*a31*a44*a52 - a15*a23*a32*a41*a54 + a15*a23*a32*a44*a51 + a15*a23*a34*a41*a52 - a15*a23*a34*a42*a51
             - a15*a24*a31*a42*a53 + a15*a24*a31*a43*a52 + a15*a24*a32*a41*a53 - a15*a24*a32*a43*a51 - a15*a24*a33*a41*a52 + a15*a24*a33*a42*a51);

    dstype c11 = (a22*a33*a44*a55 - a22*a33*a45*a54 - a22*a34*a43*a55 + a22*a34*a45*a53 + a22*a35*a43*a54 - a22*a35*a44*a53) / detA;
    dstype c21 = (a21*a33*a44*a55 - a21*a33*a45*a54 - a21*a34*a43*a55 + a21*a34*a45*a53 + a21*a35*a43*a54 - a21*a35*a44*a53) / detA;
    dstype c31 = (a21*a32*a44*a55 - a21*a32*a45*a54 - a21*a34*a42*a55 + a21*a34*a45*a52 + a21*a35*a42*a54 - a21*a35*a44*a52) / detA;
    dstype c41 = (a21*a32*a43*a55 - a21*a32*a45*a53 - a21*a33*a42*a55 + a21*a33*a45*a52 + a21*a35*a42*a53 - a21*a35*a43*a52) / detA;
    dstype c51 = (a21*a32*a43*a54 - a21*a32*a44*a53 - a21*a33*a42*a54 + a21*a33*a44*a52 + a21*a34*a42*a53 - a21*a34*a43*a52) / detA;
    dstype c12 = (a12*a33*a44*a55 - a12*a33*a45*a54 - a12*a34*a43*a55 + a12*a34*a45*a53 + a12*a35*a43*a54 - a12*a35*a44*a53) / detA;
    dstype c22 = (a11*a33*a44*a55 - a11*a33*a45*a54 - a11*a34*a43*a55 + a11*a34*a45*a53 + a11*a35*a43*a54 - a11*a35*a44*a53) / detA;
    dstype c32 = (a11*a32*a44*a55 - a11*a32*a45*a54 - a11*a34*a42*a55 + a11*a34*a45*a52 + a11*a35*a42*a54 - a11*a35*a44*a52) / detA;
    dstype c42 = (a11*a32*a43*a55 - a11*a32*a45*a53 - a11*a33*a42*a55 + a11*a33*a45*a52 + a11*a35*a42*a53 - a11*a35*a43*a52) / detA;
    dstype c52 = (a11*a32*a43*a54 - a11*a32*a44*a53 - a11*a33*a42*a54 + a11*a33*a44*a52 + a11*a34*a42*a53 - a11*a34*a43*a52) / detA;
    dstype c13 = (a12*a23*a44*a55 - a12*a23*a45*a54 - a12*a24*a43*a55 + a12*a24*a45*a53 + a12*a25*a43*a54 - a12*a25*a44*a53) / detA;
    dstype c23 = (a11*a23*a44*a55 - a11*a23*a45*a54 - a11*a24*a43*a55 + a11*a24*a45*a53 + a11*a25*a43*a54 - a11*a25*a44*a53) / detA;
    dstype c33 = (a11*a22*a44*a55 - a11*a22*a45*a54 - a11*a24*a42*a55 + a11*a24*a45*a52 + a11*a25*a42*a54 - a11*a25*a44*a52) / detA;
    dstype c43 = (a11*a22*a43*a55 - a11*a22*a45*a53 - a11*a23*a42*a55 + a11*a23*a45*a52 + a11*a25*a42*a53 - a11*a25*a43*a52) / detA;
    dstype c53 = (a11*a22*a43*a54 - a11*a22*a44*a53 - a11*a23*a42*a54 + a11*a23*a44*a52 + a11*a24*a42*a53 - a11*a24*a43*a52) / detA;
    dstype c14 = (a12*a23*a34*a55 - a12*a23*a35*a54 - a12*a24*a33*a55 + a12*a24*a35*a53 + a12*a25*a33*a54 - a12*a25*a34*a53) / detA;
    dstype c24 = (a11*a23*a34*a55 - a11*a23*a35*a54 - a11*a24*a33*a55 + a11*a24*a35*a53 + a11*a25*a33*a54 - a11*a25*a34*a53) / detA;
    dstype c34 = (a11*a22*a34*a55 - a11*a22*a35*a54 - a11*a24*a32*a55 + a11*a24*a35*a52 + a11*a25*a32*a54 - a11*a25*a34*a52) / detA;
    dstype c44 = (a11*a22*a33*a55 - a11*a22*a35*a53 - a11*a23*a32*a55 + a11*a23*a35*a52 + a11*a25*a32*a53 - a11*a25*a33*a52) / detA;
    dstype c54 = (a11*a22*a33*a54 - a11*a22*a34*a53 - a11*a23*a32*a54 + a11*a23*a34*a52 + a11*a24*a32*a53 - a11*a24*a33*a52) / detA;
    dstype c15 = (a12*a23*a34*a45 - a12*a23*a35*a44 - a12*a24*a33*a45 + a12*a24*a35*a43 + a12*a25*a33*a44 - a12*a25*a34*a43) / detA;
    dstype c25 = (a11*a23*a34*a45 - a11*a23*a35*a44 - a11*a24*a33*a45 + a11*a24*a35*a43 + a11*a25*a33*a44 - a11*a25*a34*a43) / detA;
    dstype c35 = (a11*a22*a34*a45 - a11*a22*a35*a44 - a11*a24*a32*a45 + a11*a24*a35*a42 + a11*a25*a32*a44 - a11*a25*a34*a42) / detA;
    dstype c45 = (a11*a22*a33*a45 - a11*a22*a35*a43 - a11*a23*a32*a45 + a11*a23*a35*a42 + a11*a25*a32*a43 - a11*a25*a33*a42) / detA;
    dstype c55 = (a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42) / detA;

    if (nrhs == 1) {
      dstype b1 = b[i + N*0];
      dstype b2 = b[i + N*1];
      dstype b3 = b[i + N*2];
      dstype b4 = b[i + N*3];
      dstype b5 = b[i + N*4];

      b[i + N*0] = c11*b1 + c12*b2 + c13*b3 + c14*b4 + c15*b5;
      b[i + N*1] = c21*b1 + c22*b2 + c23*b3 + c24*b4 + c25*b5;
      b[i + N*2] = c31*b1 + c32*b2 + c33*b3 + c34*b4 + c35*b5;
      b[i + N*3] = c41*b1 + c42*b2 + c43*b3 + c44*b4 + c45*b5;
      b[i + N*4] = c51*b1 + c52*b2 + c53*b3 + c54*b4 + c55*b5;        
    }
    else {
      for (int j=0; j<nrhs; j++) {
        dstype b1 = b[i + N*0 + N*5*j];
        dstype b2 = b[i + N*1 + N*5*j];     
        dstype b3 = b[i + N*2 + N*5*j];   
        dstype b4 = b[i + N*3 + N*5*j];
        dstype b5 = b[i + N*4 + N*5*j];
        b[i + N*0 + N*5*j] = c11*b1 + c12*b2 + c13*b3 + c14*b4 + c15*b5;
        b[i + N*1 + N*5*j] = c21*b1 + c22*b2 + c23*b3 + c24*b4 + c25*b5;
        b[i + N*2 + N*5*j] = c31*b1 + c32*b2 + c33*b3 + c34*b4 + c35*b5;
        b[i + N*3 + N*5*j] = c41*b1 + c42*b2 + c43*b3 + c44*b4 + c45*b5;
        b[i + N*4 + N*5*j] = c51*b1 + c52*b2 + c53*b3 + c54*b4 + c55*b5;
      }              
    }
  });
}

void GetElemNodes(dstype* unView, const dstype* uView, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
{
    int nn = np * (e2 - e1);
    int ncu = nc2 - nc1;
    int N = nn * ncu;
    int K = np * nc;

    Kokkos::parallel_for("GetElemNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx % nn;  // [0, np*(e2-e1)]
        int j = idx / nn;  // [0, ncu]
        int k = i % np;    // [0, np]
        int e = i / np + e1;
        unView[idx] = uView[k + (j + nc1) * np + e * K];
    });
}

void PutElemNodes(dstype* u, const dstype* un, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
{
    int nn = np * (e2 - e1);
    int ncu = nc2 - nc1;
    int N = nn * ncu;
    int K = np * nc;
    Kokkos::parallel_for("PutElemNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx % nn;  // [0, np*(e2-e1)]
        int j = idx / nn;  // [0, ncu]
        int k = i % np;    // [0, np]
        int e = i / np + e1;
        u[k + (j + nc1) * np + e * K] = un[idx];
    });
}

void GetFaceNodes(dstype* uh, const dstype* udg, const int* facecon, const int npf, const int ncu, const int npe, const int nc, const int f1, const int f2, const int opts) 
{
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("GetFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ndf; 
        int j = idx/ndf; // [0, ncu)
        int m = npf*f1+i;
        int k1 = facecon[2*m]; // 2*m = 2*(npf*f1+i)
        int k2 = facecon[2*m+1];
        int m1 = k1%npe; // [0, npe)
        int m2 = k2%npe; // [0, npe)
        int n1 = (k1-m1)/npe; // [0, ne)
        int n2 = (k2-m2)/npe; // [0, ne)              
        // uh npf*nf*ncu 
        if (opts==0) 
            uh[idx] = 0.5*(udg[m1+j*npe+n1*M]+udg[m2+j*npe+n2*M]);
        else if (opts==1) 
            uh[idx] = udg[m1+j*npe+n1*M];
        else if (opts==2) 
            uh[idx] = udg[m2+j*npe+n2*M];
    });
}

void GetFaceNodes(dstype* uh, const dstype* udg, const int* f2e, const int* perm, const int npf, const int ncu, const int npe, const int nc, const int nf) 
{
    int N = npf*nf*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("GetFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ncu; // [0, ncu) 
        int j = idx/ncu; // [0, npf*nf)
        int n = j%npf;   // [0, npf)
        int k = j/npf;   // [0, nf)
        int e1 = f2e[4*k+0];
        int l1 = f2e[4*k+1];
        int m = perm[n + npf*l1];
        uh[idx] = udg[m+i*npe+e1*M];
    });
}

void GetBoudaryNodes(dstype* ub, const dstype* uh, const int* boufaces, const int* elemcon, const int nfe, const int npf, const int ncu, const int nf) 
{
    int K = npf*nf;
    int N = K*ncu;
    int ndf = npf*nfe;   
    Kokkos::parallel_for("GetBoudaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx/K; // [0, ncu) 
        int j = idx%K; // [0, npf*nf)
        int n = j%npf;   // [0, npf)
        int k = j/npf;   // [0, nf)
        int e1 = boufaces[k]/nfe;
        int l1 = boufaces[k]%nfe;
        ub[idx] = uh[i + elemcon[n + npf*l1 + ndf*e1]*ncu];
    });
}

void GetBoudaryNodes(dstype* ub, const dstype* udg, const int* boufaces, const int* perm, const int nfe, const int npf, const int npe, const int ncu, const int nc, const int nf) 
{
    int K = npf*nf;
    int N = K*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("GetBoudaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx/K; // [0, ncu) 
        int j = idx%K; // [0, npf*nf)
        int n = j%npf;   // [0, npf)
        int k = j/npf;   // [0, nf)
        int e1 = boufaces[k]/nfe;
        int l1 = boufaces[k]%nfe;
        int m = perm[n + npf*l1];
        ub[idx] = udg[m+i*npe+e1*M];
    });
}

void PutBoudaryNodes(dstype* udg, const dstype* ub, const int* boufaces, const int* perm, const int nfe, const int npf, const int npe, const int ncu, const int nc, const int nf) 
{
    int K = npf*nf;
    int N = K*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("PutBoudaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx/K; // [0, ncu) 
        int j = idx%K; // [0, npf*nf)
        int n = j%npf;   // [0, npf)
        int k = j/npf;   // [0, nf)
        int e1 = boufaces[k]/nfe;
        int l1 = boufaces[k]%nfe;
        int m = perm[n + npf*l1];
        udg[m+i*npe+e1*M] = ub[idx];
    });
}

void PutBoudaryNodes(dstype* uh, const dstype* ub, const int* boufaces, const int* faceperm, const int* comperm, const int nfe, const int npf, const int ncu, const int nc, const int nf) 
{
    int K = npf*nf;
    int N = K*ncu;
    int M = nc*npf;
    int P = nc*npf*nfe;
    Kokkos::parallel_for("PutBoudaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ncu; // [0, ncu) 
        int j = idx/ncu; // [0, npf*nf)
        int n = j%npf;   // [0, npf)
        int k = j/npf;   // [0, nf)
        int e1 = boufaces[k]/nfe;
        int l1 = boufaces[k]%nfe;
        uh[comperm[i] + nc*n + M*l1 + P*e1] += ub[i + ncu*faceperm[n] + ncu*npf*k];
    });
}

void GetElementFaceNodes(dstype* uh, const dstype* udg, const int* perm, const int ndf, const int ncu, const int npe, const int nc, const int e1, const int e2) 
{
    int ne = e2-e1;
    int N = ndf*ne*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("GetElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ndf; // [0, npf*nfe] 
        int k = idx/ndf; // [0, ne*ncu)
        int e = k%ne;    // [0, ne)
        int j = k/ne;    // [0, ncu)        
        int n = perm[i]; // [0, npe)
        uh[idx] = udg[n + j*npe + (e+e1)*M];
    });
}

void GetElementFaceNodes(dstype* uhe, const dstype* uhf, const int* elemcon, const int ndf, const int ncu, const int e1, const int e2, const int opt) 
{
    int ne = e2-e1;
    int N = ndf*ne*ncu;
    if (opt==0) { // uhe = npf*nfe*ne*ncu
      Kokkos::parallel_for("GetElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
          int i = idx%ndf; // [0, npf*nfe]         
          int k = idx/ndf; // [0, ne*ncu)
          int e = k%ne;    // [0, ne)
          int j = k/ne;    // [0, ncu)                
          uhe[idx] = uhf[j + elemcon[i + (e+e1)*ndf]*ncu];
      });
    } 
    else if (opt==1) { // uhe = npf*nfe*ncu*ne
      Kokkos::parallel_for("GetElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
          int i = idx%ndf; // [0, npf*nfe]         
          int k = idx/ndf; // [0, ne*ncu)
          int j = k%ncu;    // [0, ncu)
          int e = k/ncu;    // [0, ne)                
          uhe[idx] = uhf[j + elemcon[i + (e+e1)*ndf]*ncu];
      });      
    }
    else if (opt==2) { // uhe = ncu*npf*nfe*ne
      Kokkos::parallel_for("GetElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
          int j = idx%ncu; // [0, ncu]         
          int k = idx/ncu; // [0, ndf*ne)
          int i = k%ndf;    // [0, ndf)
          int e = k/ndf;    // [0, ne)                
          uhe[idx] = uhf[j + elemcon[i + (e+e1)*ndf]*ncu];
      });      
    }
}

void PutElementFaceNodes(dstype* uhf, const dstype* uhe, const int* f2e, const int npf, const int nfe, const int ncu, const int nf) 
{
    int L = ncu*npf;
    int M = ncu*npf*nfe;
    int N = ncu*npf*nf;
    Kokkos::parallel_for("PutElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu; // [0, ncu] 
        int k = idx/ncu; // [0, npf*nf)
        int i = k%npf;    // [0, npf)
        int f = k/npf;    // [0, nf)        
        int e1 = f2e[0 + 4*f];
        int l1 = f2e[1 + 4*f];
        uhf[idx] = uhe[m + ncu*i + L*l1 + M*e1];
    });
}

void PutElementFaceNodes(dstype* uhf, const dstype* uhe, const int* f2e, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
{
    int L = ncu*npf;
    int M = ncu*npf*nfe;
    int N = ncu*npf*nf;
    Kokkos::parallel_for("PutElementFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu; // [0, ncu] 
        int k = idx/ncu; // [0, npf*nf)
        int i = k%npf;    // [0, npf)
        int f = k/npf;    // [0, nf)        
        int e2 = f2e[2 + 4*f];
        if (e2 >= 0) {
          int l2 = f2e[3 + 4*f];
          int ii = elcon[i + npf*l2 + npf*nfe*e2] - npf*f;    
          uhf[idx] += uhe[m + ncu*ii + L*l2 + M*e2];
        }
    });
}

void BlockJacobi(dstype* BE, const dstype* AE, const int* f2e, const int npf, const int nfe, const int ncu, const int nf) 
{
    int ncf = ncu*npf;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int N = ncf*ncf*nf;
    Kokkos::parallel_for("BlockJacobi", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*nf)
        int n = k%ncf;    // [0, ncf)
        int f = k/ncf;    // [0, nf)        
        int e1 = f2e[0 + 4*f];
        int l1 = f2e[1 + 4*f];
        BE[idx] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1];
    });
}

void BlockJacobi(dstype* BE, const dstype* AE, const int* f2e, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
{
    int ncf = ncu*npf;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int N = ncf*ncf*nf;
    Kokkos::parallel_for("BlockJacobi", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*nf)
        int n = k%ncf;    // [0, ncf)
        int f = k/ncf;    // [0, nf)       
        int e2 = f2e[2 + 4*f];
        if (e2 >= 0) {        
          int l2 = f2e[3 + 4*f];          
          int am = m%ncu;    // [0, ncu)
          int bm = m/ncu;   // [0, npf)
          int cm = elcon[bm + npf*l2 + npf*nfe*e2] - npf*f;    
          int dm = am + ncu*cm; // [0, ncf)
          int an = n%ncu;    // [0, ncu)          
          int bn = n/ncu;   // [0, npf)
          int cn = elcon[bn + npf*l2 + npf*nfe*e2] - npf*f;              
          int dn = an + ncu*cn; // [0, ncf)
          BE[idx] += AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2];
        }        
    });
}

void ElementalAdditiveSchwarz(dstype* BE, const dstype* AE, const int* f2e, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
{
    int ncf = ncu*npf;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int N = ncf*ncf*nf;
    Kokkos::parallel_for("ElementalAdditiveSchwarz", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*nf)
        int n = k%ncf;    // [0, ncf)
        int f = k/ncf;    // [0, nf)       
        int e2 = f2e[2 + 4*f];
        if (e2 >= 0) {        
          int e1 = f2e[0 + 4*f];
          int l1 = f2e[1 + 4*f];
          int l2 = f2e[3 + 4*f];          
          int am = m%ncu;    // [0, ncu)
          int bm = m/ncu;   // [0, npf)
          int cm = elcon[bm + npf*l2 + npf*nfe*e2] - npf*f;    
          int dm = am + ncu*cm; // [0, ncf)
          int an = n%ncu;    // [0, ncu)          
          int bn = n/ncu;   // [0, npf)
          int cn = elcon[bn + npf*l2 + npf*nfe*e2] - npf*f;              
          int dn = an + ncu*cn; // [0, ncf)
          dstype tm = AE[m + ncf*l1 + M*n + P*l1 + Q*e1] + AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2];
          BE[m + ncf*l1 + M*n + P*l1 + Q*e1] = tm;
          BE[dm + ncf*l2 + M*dn + P*l2 + Q*e2] = tm;          
        }        
    });
}

void PathAdditiveSchwarz(dstype* A, dstype* B1, dstype* C1, dstype* B2, dstype* C2, dstype* D1, dstype* D2,
        dstype* DL, dstype* DU, const dstype* AE, const int* f2e, const int* elcon, const int* epath, 
        const int* fpath, const int* lpath, const int* fintf, const int* lintf, const int npf, 
        const int nfe, const int ncu, const int ne) 
{
    int nintf = nfe - 2;
    int ncf = ncu*npf;
    int ndf = npf*nfe;
    int R = ncf*ncf;
    int K = ncf*nintf;
    int S = K*K;
    int W = K*ncf;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int N = ncf*ncf*ne;    
    Kokkos::parallel_for("PathAdditiveSchwarz", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*ne)
        int n = k%ncf;    // [0, ncf)
        int i = k/ncf;    // [0, ne)       
        int e = epath[i];
        
        int am = m%ncu;    // [0, ncu)
        int bm = m/ncu;   // [0, npf)
        int an = n%ncu;    // [0, ncu)          
        int bn = n/ncu;   // [0, npf)        
        int e1, l1, f1, e2, l2, f2, cm, dm, cn, dn;
        for (int i1=0; i1<nintf; i1++)
        {
          f1 = fintf[i1 + nintf*i];
          e1 = f2e[0 + 4*f1];
          l1 = f2e[1 + 4*f1];
          e2 = f2e[2 + 4*f1];
          if (e2>=0) {
            l2 = f2e[3 + 4*f1];      
            cm = elcon[bm + npf*l2 + ndf*e2] - npf*f1;    
            dm = am + ncu*cm; // [0, ncf)            
            cn = elcon[bn + npf*l2 + ndf*e2] - npf*f1;              
            dn = an + ncu*cn; // [0, ncf)
            A[m + ncf*i1 + K*n + W*i1 + S*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1] + AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2]; 
          }
          else {
            A[m + ncf*i1 + K*n + W*i1 + S*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1];
          }
          
          l1 = lintf[i1 + nintf*i];
          cm = elcon[bm + npf*l1 + ndf*e] - npf*f1;    
          dm = am + ncu*cm; // [0, ncf)
          for (int i2=0; i2<nintf; i2++) {
            if (i1 != i2) {
              l2 = lintf[i2 + nintf*i];
              f2 = fintf[i2 + nintf*i];
              cn = elcon[bn + npf*l2 + ndf*e] - npf*f2;              
              dn = an + ncu*cn; // [0, ncf)              
              A[m + ncf*i1 + K*n + W*i2 + S*i] = AE[dm + ncf*l1 + M*dn + P*l2 + Q*e];
            }
          }
          
          l2 = lpath[0 + 2*i];
          f2 = fpath[0 + 2*i];
          cn = elcon[bn + npf*l2 + ndf*e] - npf*f2;              
          dn = an + ncu*cn; // [0, ncf)                        
          B1[m + ncf*i1 + K*n + W*i] = AE[dm + ncf*l1 + M*dn + P*l2 + Q*e];
          C1[m + ncf*n + R*i1 + W*i] = AE[dn + ncf*l2 + M*dm + P*l1 + Q*e];
          
          l2 = lpath[1 + 2*i];
          f2 = fpath[1 + 2*i];
          cn = elcon[bn + npf*l2 + ndf*e] - npf*f2;              
          dn = an + ncu*cn; // [0, ncf)                        
          B2[m + ncf*i1 + K*n + W*i] = AE[dm + ncf*l1 + M*dn + P*l2 + Q*e];
          C2[m + ncf*n + R*i1 + W*i] = AE[dn + ncf*l2 + M*dm + P*l1 + Q*e];    
        }
        
        l1 = lpath[0 + 2*i];
        f1 = fpath[0 + 2*i];
        cm = elcon[bm + npf*l1 + npf*nfe*e] - npf*f1;    
        dm = am + ncu*cm; // [0, ncf)
        
        l2 = lpath[1 + 2*i];
        f2 = fpath[1 + 2*i];
        cn = elcon[bn + npf*l2 + npf*nfe*e] - npf*f2;              
        dn = an + ncu*cn; // [0, ncf)                        
                
        DU[m + ncf*n + R*i] = AE[dm + ncf*l1 + M*dn + P*l2 + Q*e];
        DL[m + ncf*n + R*i] = AE[dn + ncf*l2 + M*dm + P*l1 + Q*e];
        
        e1 = f2e[0 + 4*f1];
        l1 = f2e[1 + 4*f1];
        e2 = f2e[2 + 4*f1];
        if (e2>=0) {
          l2 = f2e[3 + 4*f1];      
          cm = elcon[bm + npf*l2 + ndf*e2] - npf*f1;    
          dm = am + ncu*cm; // [0, ncf)            
          cn = elcon[bn + npf*l2 + ndf*e2] - npf*f1;              
          dn = an + ncu*cn; // [0, ncf)            
          D1[m + ncf*n + R*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1] + AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2]; 
        }
        else {
          D1[m + ncf*n + R*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1];
        }
        
        e1 = f2e[0 + 4*f2];
        l1 = f2e[1 + 4*f2];
        e2 = f2e[2 + 4*f2];
        if (e2>=0) {
          l2 = f2e[3 + 4*f2];      
          cm = elcon[bm + npf*l2 + ndf*e2] - npf*f1;    
          dm = am + ncu*cm; // [0, ncf)            
          cn = elcon[bn + npf*l2 + ndf*e2] - npf*f1;              
          dn = an + ncu*cn; // [0, ncf)                    
          D2[m + ncf*n + R*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1] + AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2]; 
        }
        else {
          D2[m + ncf*n + R*i] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1];
        }
    });
}

void AssembleBlockILU0(dstype* BE, const dstype* AE, const int* f2e, const int* elcon, const int* face, const int* row_ptr, const int* col_ind, const int npf, const int nfe, const int ncu, const int nf, const int nb) 
{
    int ncf = ncu*npf;
    int ndf = npf*nfe;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int R = ncf*ncf;
    int S = R*nb;
    int N = S*nf;
    Kokkos::parallel_for("AssembleCRSMatrix", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*nb*nf)
        int n = k%ncf;   // [0, ncf)
        int q = k/ncf;   // [0, nb*nf)        
        int r = q%nb;    // [0, nb)        
        int i = q/nb;    // [0, nf)        
        
        int am = m%ncu;    // [0, ncu)
        int bm = m/ncu;   // [0, npf)
        int an = n%ncu;    // [0, ncu)          
        int bn = n/ncu;   // [0, npf)        
        
        int row_start = row_ptr[i];
        int row_end = row_ptr[i+1];    
        
        int fi = face[r + nb*i];
        int e1 = f2e[0 + 4*fi];
        int l1 = f2e[1 + 4*fi];
        int e2 = f2e[2 + 4*fi];
        int l2 = f2e[3 + 4*fi];
        
        int m1, n1, m2, n2;
        int nfi = npf*fi;
        m1 = am + ncu*(elcon[bm + npf*l1 + ndf*e1] - nfi);
        n1 = an + ncu*(elcon[bn + npf*l1 + ndf*e1] - nfi);                
        if (e2>=0) {
          m2 = am + ncu*(elcon[bm + npf*l2 + ndf*e2] - nfi);
          n2 = an + ncu*(elcon[bn + npf*l2 + ndf*e2] - nfi);                
          BE[m + ncf*n + R*r + S*row_start] = AE[m1 + ncf*l1 + M*n1 + P*l1 + Q*e1] + AE[m2 + ncf*l2 + M*n2 + P*l2 + Q*e2]; 
        }
        else {
          BE[m + ncf*n + R*r + S*row_start] = AE[m1 + ncf*l1 + M*n1 + P*l1 + Q*e1];
        }
        
        for(int t = row_start+1; t<row_end; t++)
        {
          int j = col_ind[t];  
          int fj = face[r + nb*j];
          int je1 = f2e[0 + 4*fj];       
          int je2 = f2e[2 + 4*fj];     
          int e=0, k1=0, k2=0;
          if (je1 == e1) {
            e = e1;                    
            k1 = l1;
            k2 = f2e[1 + 4*fj];
          }
          else if (je1 == e2) {
            e = e2;                    
            k1 = l2;
            k2 = f2e[1 + 4*fj];    
          }
          else if (je2 == e1) {
            e = e1;                    
            k1 = l1;
            k2 = f2e[3 + 4*fj];     
          }
          else if (je2 == e2) {
            e = e2;                    
            k1 = l2;
            k2 = f2e[3 + 4*fj];
          }          
          m1 = am + ncu*(elcon[bm + npf*k1 + ndf*e] - nfi);
          n2 = an + ncu*(elcon[bn + npf*k2 + ndf*e] - npf*fj);
          BE[m + ncf*n + R*r + S*t] = AE[m1 + ncf*k1 + M*n2 + P*k2 + Q*e];
        }                        
    });
}

// void AssembleResidual(dstype* R, const dstype* Rh, const int* e2f, const int* elcon, const int npf, const int nfe, const int ncu, const int ne) 
// {
//     //int nfe2 = 2*(nfe-1);
//     int ncf = ncu*npf;
//     int M = ncf*nfe;
//     int N = M*ne;
//     Kokkos::parallel_for("AssembleResidual", N, KOKKOS_LAMBDA(const size_t idx) {
//         int m = idx%ncf; // [0, ncf] 
//         int q = idx/ncf; // [0, nfe*ne)
//         int k = q%nfe;   // [0, nfe)
//         int e = q/nfe;   // [0, ne)                
//         int fk = e2f[k + nfe*e];
//         int am = m%ncu;    // [0, ncu)
//         int bm = m/ncu;    // [0, npf)
//         int cm = elcon[bm + npf*k + npf*nfe*e] - npf*fk;    
//         int dm = am + ncu*cm; // [0, ncf)            
//         Kokkos::atomic_add(&R[m + ncf*fk], Rh[dm + ncf*k + M*e]);
//     });
// }

// void AssembleJacobianMatrix(dstype* BE, const dstype* AE, const int* e2f, const int* f2f, const int* elcon, const int npf, const int nfe, const int ncu, const int ne) 
// {
//     int nfe2 = 2*(nfe-1);
//     int ncf = ncu*npf;
//     int M = ncf*nfe;
//     int P = M*ncf;
//     int Q = M*M;
//     int R = ncf*ncf;
//     int N = Q*ne;
//     Kokkos::parallel_for("AssembleJacobianMatrix", N, KOKKOS_LAMBDA(const size_t idx) {
//         int m = idx%ncf; // [0, ncf] 
//         int q = idx/ncf; // [0, nfe*ncf*nfe*ne)
//         int k = q%nfe;   // [0, nfe)
//         int s = q/nfe;   // [0, ncf*nfe*ne)                
//         int n = s%ncf;   // [0, ncf)
//         int t = s/ncf;   // [0, nfe*ne)        
//         int l = t%nfe;   // [0, nfe)        
//         int e = t/nfe;   // [0, ne)                
//         int fk = e2f[k + nfe*e];
//         int fl = e2f[l + nfe*e];
//         int j = 0;                
//         if (k != l) {
//           for (int i=0; i<nfe2; i++)
//             if (f2f[i + nfe2*fk] == fl)
//               j = i + 1;
//         }       
//         int am = m%ncu;    // [0, ncu)
//         int bm = m/ncu;    // [0, npf)
//         int cm = elcon[bm + npf*k + npf*nfe*e] - npf*fk;    
//         int dm = am + ncu*cm; // [0, ncf)            
//         int an = n%ncu;    // [0, ncu)          
//         int bn = n/ncu;    // [0, npf)
//         int cn = elcon[bn + npf*l + npf*nfe*e] - npf*fl;      
//         int dn = an + ncu*cn; // [0, ncf)       
//         if (j==0) 
//           Kokkos::atomic_add(&BE[m + ncf*n + R*(j + (2*nfe-1)*fk)], AE[dm + ncf*k + M*dn + P*l + Q*e]);
//         else
//           BE[m + ncf*n + R*(j + (2*nfe-1)*fk)] = AE[dm + ncf*k + M*dn + P*l + Q*e];                
//     });
// }

void AssembleJacobianMatrix(dstype* BE, const dstype* AE, const int* f2e, const int* f2f, const int* f2l, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
{
    int nfe2 = 2*(nfe-1);
    int ncf = ncu*npf;
    int M = ncf*nfe;
    int P = M*ncf;
    int Q = M*M;
    int R = ncf*ncf;
    int N = R*nf;
    Kokkos::parallel_for("AssembleJacobianMatrix", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, ncf*nf)
        int n = k%ncf;    // [0, ncf)
        int f = k/ncf;    // [0, nf)        
        int e1 = f2e[0 + 4*f];
        int l1 = f2e[1 + 4*f];
        int e2 = f2e[2 + 4*f];
        int l2 = f2e[3 + 4*f];
        int an = n%ncu;    // [0, ncu)          
        int bn = n/ncu;   // [0, npf)        
        BE[m + ncf*n + R*((2*nfe-1)*f)] = AE[m + ncf*l1 + M*n + P*l1 + Q*e1]; 
        for (int j=0; j<nfe-1; j++) {
          int lj = f2l[j + nfe2*f];          
          int fj = f2f[j + nfe2*f];          
          int cn = elcon[bn + npf*lj + npf*nfe*e1] - npf*fj;      
          int dn = an + ncu*cn; // [0, ncf)
          BE[m + ncf*n + R*(1 + j + (2*nfe-1)*f)] = AE[m + ncf*l1 + M*dn + P*lj + Q*e1];
        }
        if (e2 >= 0) {  
          int am = m%ncu;    // [0, ncu)
          int bm = m/ncu;   // [0, npf)
          int cm = elcon[bm + npf*l2 + npf*nfe*e2] - npf*f;    
          int dm = am + ncu*cm; // [0, ncf)          
          int cn = elcon[bn + npf*l2 + npf*nfe*e2] - npf*f;      
          int dn = an + ncu*cn; // [0, ncf)
          BE[m + ncf*n + R*((2*nfe-1)*f)] += AE[dm + ncf*l2 + M*dn + P*l2 + Q*e2];
          for (int j=0; j<nfe-1; j++) {            
            int lj = f2l[nfe-1 + j + nfe2*f];
            int fj = f2f[nfe-1 + j + nfe2*f];
            cn = elcon[bn + npf*lj + npf*nfe*e2] - npf*fj;      
            dn = an + ncu*cn; // [0, ncf)
            BE[m + ncf*n + R*(nfe + j + (2*nfe-1)*f)] = AE[dm + ncf*l2 + M*dn + P*lj + Q*e2];
          }
        }        
    });
}

// void BlockJacobian1(dstype* BE, const dstype* AE, const int* f2e, const int* f2f, const int* f2l, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
// {
//     int nfe2 = 2*(nfe-1);
//     int ncf = ncu*npf;
//     int M = ncf*nfe;
//     int P = M*ncf;
//     int Q = M*M;
//     int R = ncf*ncf;
//     int N = R*nf;
//     Kokkos::parallel_for("BlockJacobian1", N, KOKKOS_LAMBDA(const size_t idx) {
//         int m = idx%ncf; // [0, ncf] 
//         int k = idx/ncf; // [0, ncf*nf)
//         int n = k%ncf;    // [0, ncf)
//         int f = k/ncf;    // [0, nf)        
//         int e1 = f2e[0 + 4*f];
//         int l1 = f2e[1 + 4*f];
//         for (int j=0; j<nfe-1; j++) {
//           int lj = f2l[j + nfe2*f];          
//           int fj = f2f[j + nfe2*f];          
//           int an = n%ncu;    // [0, ncu)          
//           int bn = n/ncu;   // [0, npf)
//           int cn = elcon[bn + npf*lj + npf*nfe*e1] - npf*fj;      
//           int dn = an + ncu*cn; // [0, ncf)
//           BE[m + ncf*n + R*(j + (nfe-1)*f)] = AE[m + ncf*l1 + M*dn + P*lj + Q*e1];
//         }
//     });
// }
// 
// void BlockJacobian2(dstype* BE, const dstype* AE, const int* f2e, const int* f2f, const int* f2l, const int* elcon, const int npf, const int nfe, const int ncu, const int nf) 
// {
//     int nfe2 = 2*(nfe-1);
//     int ncf = ncu*npf;
//     int M = ncf*nfe;
//     int P = M*ncf;
//     int Q = M*M;
//     int R = ncf*ncf;
//     int N = R*nf;
//     Kokkos::parallel_for("BlockJacobian1", N, KOKKOS_LAMBDA(const size_t idx) {
//         int m = idx%ncf; // [0, ncf] 
//         int k = idx/ncf; // [0, ncf*nf)
//         int n = k%ncf;    // [0, ncf)
//         int f = k/ncf;    // [0, nf)        
//         int e2 = f2e[2 + 4*f];
//         int l2 = f2e[3 + 4*f];
//         if (e2 >= 0) {  
//           for (int j=0; j<nfe-1; j++) {
//             int am = m%ncu;    // [0, ncu)
//             int bm = m/ncu;   // [0, npf)
//             int cm = elcon[bm + npf*l2 + npf*nfe*e2] - npf*f;    
//             int dm = am + ncu*cm; // [0, ncf)
//             int lj = f2l[nfe-1 + j + nfe2*f];
//             int fj = f2f[nfe-1 + j + nfe2*f];
//             int an = n%ncu;    // [0, ncu)          
//             int bn = n/ncu;   // [0, npf)
//             int cn = elcon[bn + npf*lj + npf*nfe*e2] - npf*fj;      
//             int dn = an + ncu*cn; // [0, ncf)
//             BE[m + ncf*n + R*(j + (nfe-1)*f)] = AE[dm + ncf*l2 + M*dn + P*lj + Q*e2];
//             //printf("%d %d %d %d %d %d %d %g\n", dm, l2, dn, cn, fj, lj, e2, AE[dm + ncf*l2 + M*dn + P*lj + Q*e2]);
//           }
//         }
//     });
// }

void ApplyFace2Face(dstype* Rf, const dstype* Rh, const int* f2f, const int npf, const int nfe, const int ncu, const int nf, const int offset) 
{
    int nfe1 = nfe-1;
    int nfe2 = 2*(nfe-1);
    int ncf = ncu*npf;    
    int N = ncf*(nfe-1)*nf;     
    Kokkos::parallel_for("ApplyFace2Face", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, (nfe-1)*nf)
        int j = k%nfe1;    // [0, nfe-1)
        int f = k/nfe1;    // [0, nf)       
        int fj = f2f[offset + j + nfe2*f];       
        Rf[idx] = Rh[m + ncf*fj];                
    });
}

void ApplyFace2Face(dstype* Rf, const dstype* Rh, const int* f2f, const int npf, const int nfe, const int ncu, const int nf) 
{
    int nfe1 = 2*nfe-1;
    int nfe2 = 2*(nfe-1);
    int ncf = ncu*npf;    
    int N = ncf*(2*nfe-1)*nf;     
    Kokkos::parallel_for("ApplyFace2Face", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncf; // [0, ncf] 
        int k = idx/ncf; // [0, (2*nfe-1)*nf)
        int j = k%nfe1;    // [0, 2*nfe-1)
        int f = k/nfe1;    // [0, nf)       
        int fj = (j==0) ? f : f2f[j-1 + nfe2*f];       
        Rf[idx] = Rh[m + ncf*fj];                
    });
}

void GetBoundaryNodes(dstype* uh, const dstype* udg, const int* boufaces, const int ngf, const int nfe, const int ne, const int nc, const int nfaces) 
{
    int N = ngf*nfaces*nc;
    int M = ngf*nfe*ne;
    Kokkos::parallel_for("GetBoundaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ngf; // [0, ngf] 
        int k = idx/ngf; // [0, nfaces*nc)
        int e = k%nfaces;    // [0, nfaces)
        int j = k/nfaces;    // [0, nc)        
        int n = boufaces[e]; // [0, nfe*ne)
        uh[idx] = udg[i + n*ngf + j*M];
    });
}


void PutBoundaryNodes(dstype* udg, const dstype* uh, const int* boufaces, const int ngf, const int nfe, const int ne, const int nc, const int nfaces) 
{
    int N = ngf*nfaces*nc;
    int M = ngf*nfe*ne;
    Kokkos::parallel_for("PutBoundaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%ngf; // [0, ngf] 
        int k = idx/ngf; // [0, nfaces*nc)
        int e = k%nfaces;    // [0, nfaces)
        int j = k/nfaces;    // [0, nc)        
        int n = boufaces[e]; // [0, nfe*ne)
        udg[i + n*ngf + j*M] = uh[idx];
    });
}

void PutFaceNodes(dstype* udg, const dstype* uh, const int* facecon, const int npf, const int ncu, const int npe, const int nc, const int f1, const int f2)
{
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    int M = npe*nc;
    Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int p = idx%npf; // [0, npf)
        int q = idx/npf; // [0, ncu*nf)        
        int j = q%ncu;  // [0, ncu]
        int f = q/ncu; // [0, nf)
        int i = p + npf*f;
        int m = npf*f1+i;
        int k1 = facecon[2*m]; // 2*m = 2*(npf*f1 + p + npf*f) = 2*p + 2*npf*(f1 + f)
        int k2 = facecon[2*m+1];
        int m1 = k1%npe; // [0, npe)
        int m2 = k2%npe; // [0, npe)
        int n1 = (k1-m1)/npe; // [0, ne)
        int n2 = (k2-m2)/npe; // [0, ne)                          
        dstype uhx = uh[idx]; // npf*ncu*nf               
        if (k1==k2) { // boundary face 
            Kokkos::atomic_sub(&udg[m1+j*npe+n1*M], uhx);            
        }
        else { // interior face
            Kokkos::atomic_sub(&udg[m1+j*npe+n1*M], uhx);
            Kokkos::atomic_add(&udg[m2+j*npe+n2*M], uhx);
        }        
    });                        
}

void assembleMatrixE(dstype* E, const dstype* Etmp, const int* facecon, const int* f2e, const int npf, const int npe, const int nfe, const int f1, const int f2)
{
    int nf = f2-f1;    
    int N = npf*npf*nf;    
    Kokkos::parallel_for("assembleMatrixE", N, KOKKOS_LAMBDA(const size_t idx) {
        int k = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nf)        
        int i = q%npf;  // [0, npf]
        int j = q/npf; // [0, nf)
        int e1 = f2e[4*j+0];
        int l1 = f2e[4*j+1];
        int e2 = f2e[4*j+2];
        int l2 = f2e[4*j+3];
        int n1 = facecon[0 + 2*k + 2*npf*j];
        int n2 = facecon[1 + 2*k + 2*npf*j];
        int m1 = n1%npe; // [0, npe)
        int m2 = n2%npe; // [0, npe)
        dstype a = Etmp[k + npf*i + npf*npf*j];
        // res.E = npe*npf*nfe*ne*nd
        if (e1==e2) {
            Kokkos::atomic_sub(&E[m1 + npe*i + npe*npf*l1 + npe*npf*nfe*e1], a);
        }
        else {
            Kokkos::atomic_sub(&E[m1 + npe*i + npe*npf*l1 + npe*npf*nfe*e1], a);
            Kokkos::atomic_add(&E[m2 + npe*i + npe*npf*l2 + npe*npf*nfe*e2], a);          
        }                
    });                        
}

void assembleMatrixE(dstype* E, const dstype* Etmp, const int* perm, const int npf, const int npe, const int nfe, const int ne)
{    
    int M = npe*npf*nfe;
    int N = npf*npf*nfe*ne;    
    Kokkos::parallel_for("assembleMatrixE", N, KOKKOS_LAMBDA(const size_t idx) {
        int k = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nfe*ne)        
        int i = q%npf;  // [0, npf]
        int j = q/npf; //   [0, nfe*ne)        
        int l = j%nfe; //   [0, nfe]
        int e = j/nfe; // [0, ne]
        int m = perm[k + npf*l];
        Kokkos::atomic_add(&E[m + npe*i + npe*npf*l + M*(e)], Etmp[idx]);        
    });                        
}

void assembleRu(dstype* Ru, const dstype* Rutmp, const int* perm, const int npe, const int ndf, const int ne)
{        
    int N = ndf*ne;    
    Kokkos::parallel_for("assembleRu", N, KOKKOS_LAMBDA(const size_t idx) {
        int k = idx%ndf; // [0, ndf)                  
        int e = idx/ndf; // [0, ne]        
        int m = perm[k];
        Kokkos::atomic_sub(&Ru[m + npe*e], Rutmp[idx]);        
    });                        
}

void assembleMatrixBD(dstype* D, const dstype* Dtmp, const int* perm, const int npe, const int npf, const int nfe, const int ne)
{        
    int M = npe*npe;    
    int N = npf*npf*nfe*ne;    
    Kokkos::parallel_for("assembleMatrixBD", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nfe*ne)                  
        int j = q%npf;   // [0, npf]
        int p = q/npf;   // [0, nfe*ne)        
        int l = p%nfe;   // [0, nfe]
        int e = p/nfe;   // [0, ne]                
        int m = perm[i + npf*l];
        int n = perm[j + npf*l];
        Kokkos::atomic_add(&D[m + npe*n + M*e], Dtmp[idx]);        
    });                        
}

// npf*npf*nfe*ne*ncu*ncu -> npf*nfe*npe*ne*ncu*ncu
// assembleMatrixGK(res.K, Dtmp, mesh.perm, npe, npf, nfe, ne*ncu*ncu);
void assembleMatrixGK(dstype* K, const dstype* Ktmp, const int* perm, const int npe, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe; 
    int M = npf*nfe*npe;    
    int N = npf*npf*nfe*ne;    
    Kokkos::parallel_for("assembleMatrixGK", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nfe*ne)                  
        int j = q%npf;   // [0, npf]
        int p = q/npf;   // [0, nfe*ne)        
        int l = p%nfe;   // [0, nfe]
        int e = p/nfe;   // [0, ne]                        
        int n = perm[j + npf*l];
        Kokkos::atomic_add(&K[i + npf*l + ndf*n + M*e], Ktmp[idx]);        
    });                        
}

void assembleMatrixF(dstype* F, const dstype* Ftmp, const int* perm, const int npe, const int npf, const int nfe, const int ne)
{        
    int M = npe*npf;
    int K = M*nfe;    
    int N = npf*npf*nfe*ne;    
    Kokkos::parallel_for("assembleMatrixF", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nfe*ne)                  
        int j = q%npf;   // [0, npf]
        int p = q/npf;   // [0, nfe*ne)        
        int l = p%nfe;   // [0, nfe]
        int e = p/nfe;   // [0, ne]                
        int m = perm[i + npf*l];        
        Kokkos::atomic_add(&F[m + npe*j + M*l + K*e], Ftmp[idx]);        
    });                        
}

void assembleMatrixH(dstype* H, const dstype* Htmp, const int* perm, const int npe, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int L = ndf*npf;  
    int M = ndf*ndf;    
    int N = npf*npf*nfe*ne;    
    Kokkos::parallel_for("assembleMatrixH", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf)
        int q = idx/npf; // [0, npf*nfe*ne)                  
        int j = q%npf;   // [0, npf]
        int p = q/npf;   // [0, nfe*ne)        
        int l = p%nfe;   // [0, nfe]
        int e = p/nfe;   // [0, ne]                        
        Kokkos::atomic_add(&H[i + npf*l + ndf*j + L*l + M*e], Htmp[idx]);        
    });                        
}

// npf*npf*nfaces*ncu12*ncu -> ncu12*npf*npe*ncu*nfaces        
void assembleMatrixKint(dstype* udg, const dstype* uh, const int* boufaces, const int* perm, const int npe, const int npf, const int nfe, const int ncu12, const int ncu, const int nfaces)
{
    int N = npf*npf*nfaces*ncu12*ncu; 
    int M2 = ncu12*npf;
    int M3 = ncu12*npf*npe;
    int M4 = ncu12*npf*npe*ncu;
    Kokkos::parallel_for("PutBoundaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf] 
        int k = idx/npf; // [0, npf*nfaces*ncu12*ncu]
        int j = k%npf;   // [0, npf]
        int m = k/npf;   // [0, nfaces*ncu12*ncu]
        int e = m%nfaces;    // [0, nfaces]
        int c = m/nfaces;    // [0, ncu12*ncu]
        int a = c%ncu12; // [0, ncu12]
        int b = c/ncu12; // [0, ncu]        
        int l = boufaces[e]%nfe; // [0, nfe]        
        int n = perm[j + npf*l]; // [0, npe]
        udg[a + ncu12*i + M2*n + M3*b + M4*e] = uh[idx];
    });
}

// npf*npf*nfaces*ncu12*ncq ->  npf*npe*nfaces*ncu12*ncq
void assembleMatrixGint(dstype* udg, const dstype* uh, const int* boufaces, const int* perm, const int npe, const int npf, const int nfe, const int ncu12, const int ncq, const int nfaces)
{
    int N = npf*npf*nfaces*ncu12*ncq; 
    int M2 = npf*npe;
    int M3 = npf*npe*nfaces;
    int M4 = npf*npe*nfaces*ncu12;
    Kokkos::parallel_for("PutBoundaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf] 
        int k = idx/npf; // [0, npf*nfaces*ncu12*ncq]
        int j = k%npf;   // [0, npf]
        int m = k/npf;   // [0, nfaces*ncu12*ncq]
        int e = m%nfaces;    // [0, nfaces]
        int c = m/nfaces;    // [0, ncu12*ncq]
        int a = c%ncu12; // [0, ncu12]
        int b = c/ncu12; // [0, ncq]        
        int l = boufaces[e]%nfe; // [0, nfe]        
        int n = perm[j + npf*l]; // [0, npe]
        udg[i + npf*n + M2*e + M3*a + M4*b] = uh[idx];
    });
}

// npf*npf*nfaces*ncu12*ncu -> ncu12*npf*ncu*npf*nfe*nfaces        
void assembleMatrixHint(dstype* udg, const dstype* uh, const int* boufaces, const int npe, const int npf, const int nfe, const int ncu12, const int ncu, const int nfaces)
{
    int N = npf*npf*nfaces*ncu12*ncu; 
    int M2 = ncu12*npf;
    int M3 = ncu12*npf*ncu;
    int M4 = ncu12*npf*ncu*npf*nfe;
    Kokkos::parallel_for("PutBoundaryNodes", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npf; // [0, npf] 
        int k = idx/npf; // [0, npf*nfaces*ncu12*ncu]
        int j = k%npf;   // [0, npf]
        int m = k/npf;   // [0, nfaces*ncu12*ncu]
        int e = m%nfaces;    // [0, nfaces]
        int c = m/nfaces;    // [0, ncu12*ncu]
        int a = c%ncu12; // [0, ncu12]
        int b = c/ncu12; // [0, ncu]        
        int l = boufaces[e]%nfe; // [0, nfe]                
        udg[a + ncu12*i + M2*b + M3*(j + npf*l) + M4*e] = uh[idx];
    });
}

// npe*npe*ne*ncu*ncu -> npe*ncu*npe*ncu*ne
void schurMatrixD(dstype* D, const dstype* Dtmp,  const int npe, const int ncu, const int ne)
{        
    int M = npe*npe;    
    int L = npe*npe*ne;
    int K = npe*npe*ne*ncu;
    int N = npe*ncu*npe*ncu*ne;    
    Kokkos::parallel_for("schurMatrixD", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npe; // [0, npe)
        int q = idx/npe; // [0, ncu*npe*ncu*ne)                  
        int m = q%ncu;   // [0, ncu]
        int p = q/ncu;   // [0, npe*ncu*ne)        
        int j = p%npe;   // [0, npe]
        int s = p/npe;   // [0, ncu*ne]                
        int n = s%ncu;   // [0, ncu]
        int e = s/ncu;   // [0, ne]
        D[idx] = Dtmp[i + npe*j + M*e + L*m + K*n];
    });                        
}

// [npe*npe*ne*ncu*ncu] x [npe*npe*ne] -> npe*ncu*npe*ncu*ne
void schurMatrixBMinvC(dstype* D, const dstype *B, const dstype *MinvC, dstype scalar, const int npe, const int ncu, const int ne)
{        
    int M = npe*npe;    
    int L = npe*npe*ne;
    int K = npe*npe*ne*ncu;
    int N = npe*ncu*npe*ncu*ne;    
    Kokkos::parallel_for("schurMatrixBMinvC", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npe; // [0, npe)
        int q = idx/npe; // [0, ncu*npe*ncu*ne)                  
        int m = q%ncu;   // [0, ncu]
        int p = q/ncu;   // [0, npe*ncu*ne)        
        int j = p%npe;   // [0, npe]
        int s = p/npe;   // [0, ncu*ne]                
        int n = s%ncu;   // [0, ncu]
        int e = s/ncu;   // [0, ne]
        dstype sum = 0;
        for (int k=0; k<npe; k++)
            sum += scalar*B[i + npe*k + M*e + L*m + K*n]*MinvC[k + npe*j + M*e];
        D[idx] += sum;
    });                        
}

// npe*(npf*nfe)*ne*(ncu*ncu) -> npe*(ncu*ncu)*(npf*nfe)*ne
void schurMatrixF(dstype* F, const dstype* Ftmp,  const int npe, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ncu2 = ncu*ncu;
    int ndf = npf*nfe;
    int M = npe*ndf;    
    int L = npe*ndf*ne;    
    int N = npe*(ncu*ncu)*(npf*nfe)*ne;    
    Kokkos::parallel_for("schurMatrixF", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npe; // [0, npe)
        int q = idx/npe; // [0, (ncu*ncu)*(npf*nfe)*ne)                  
        int m = q%ncu2;   // [0, ncu*ncu]
        int p = q/ncu2;   // [0, (npf*nfe)*ne)        
        int j = p%ndf;   // [0, npf*nfe]
        int e = p/ndf;   // [0, ne]                        
        F[idx] = Ftmp[i + npe*j + M*e + L*m];
    });                        
}

// [npe*npe*ne*(ncu*ncu)] x [npe*(npf*nfe)*ne] -> npe*(ncu*ncu)*(npf*nfe)*ne
void schurMatrixBMinvE(dstype* F, const dstype *B, const dstype *MinvE, dstype scalar, const int npe, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ncu2 = ncu*ncu;
    int ndf = npf*nfe;
    int M = npe*ndf;    
    int L = npe*npe*ne;    
    int N = npe*(ncu*ncu)*(npf*nfe)*ne;      
    Kokkos::parallel_for("schurMatrixBMinvE", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npe; // [0, npe)
        int q = idx/npe; // [0, (ncu*ncu)*(npf*nfe)*ne)                  
        int m = q%ncu2;   // [0, ncu*ncu]
        int p = q/ncu2;   // [0, (npf*nfe)*ne)        
        int j = p%ndf;   // [0, npf*nfe]
        int e = p/ndf;   // [0, ne]                        
        dstype sum = 0;
        for (int k=0; k<npe; k++)
            sum += scalar*B[i + npe*k + npe*npe*e + L*m]*MinvE[k + npe*j + M*e];
        F[idx] -= sum;
    });                        
}

// (npf*nfe)*npe*ne*ncu*ncu -> ncu*(npf*nfe)*npe*ncu*ne
void schurMatrixK(dstype* K, const dstype* Ktmp,  const int npe, const int ncu12, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int M = npe*ndf;    
    int L = npe*ndf*ne; 
    int P = L*ncu12;   
    int N = ncu12*(npf*nfe)*npe*ncu*ne;    
    Kokkos::parallel_for("schurMatrixK", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu12; // [0, ncu12)
        int q = idx/ncu12; // [0, (npf*nfe)*npe*ncu12*ne)                  
        int i = q%ndf;   // [0, npf*nfe]
        int p = q/ndf;   // [0, npe*ncu*ne)        
        int j = p%npe;   // [0, npe]
        int s = p/npe;   // [0, ncu*ne]                        
        int n = s%ncu;   // [0, ncu]
        int e = s/ncu;   // [0, ne]
        K[idx] = Ktmp[i + ndf*j + M*e + L*m + P*n];
    });                        
}

// [(npf*nfe)*npe*ne*ncu*ncu] x [npe*npe*ne] -> ncu*(npf*nfe)*npe*ncu*ne
void schurMatrixGMinvC(dstype* K, const dstype *G, const dstype *MinvC, dstype scalar, const int npe, const int ncu12, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int Q = npe*npe;
    int M = npe*ndf;    
    int L = npe*ndf*ne; 
    int P = L*ncu12;   
    int N = ncu12*(npf*nfe)*npe*ncu*ne;
    Kokkos::parallel_for("schurMatrixGMinvC", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu12; // [0, ncu12)
        int q = idx/ncu12; // [0, (npf*nfe)*npe*ncu12*ne)                  
        int i = q%ndf;   // [0, npf*nfe]
        int p = q/ndf;   // [0, npe*ncu*ne)        
        int j = p%npe;   // [0, npe]
        int s = p/npe;   // [0, ncu*ne]                        
        int n = s%ncu;   // [0, ncu]
        int e = s/ncu;   // [0, ne]
        dstype sum = 0;
        for (int k=0; k<npe; k++)
            sum += scalar*G[i + ndf*k + M*e + L*m + P*n]*MinvC[k + npe*j + Q*e];
        K[idx] += sum;
    });                        
}

// (npf*nfe)*(npf*nfe)*ne*ncu*ncu -> ncu*(npf*nfe)*ncu*(npf*nfe)*ne
void schurMatrixH(dstype* H, const dstype* Htmp,  const int ncu12, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int M = ndf*ndf;    
    int L = ndf*ndf*ne; 
    int P = L*ncu12;   
    int N = ncu12*(npf*nfe)*ncu*(npf*nfe)*ne;    
    Kokkos::parallel_for("schurMatrixH", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu12; // [0, ncu12)
        int q = idx/ncu12; // [0, (npf*nfe)*npe*ncu12*ne)                  
        int i = q%ndf;   // [0, npf*nfe]
        int p = q/ndf;   // [0, ncu*(npf*nfe)*ne)        
        int n = p%ncu;   // [0, ncu]
        int s = p/ncu;   // [0, (npf*nfe)*ne]                        
        int j = s%ndf;   // [0, (npf*nfe)]
        int e = s/ndf;   // [0, ne]
        H[idx] = Htmp[i + ndf*j + M*e + L*m + P*n];
    });                        
}

// [(npf*nfe)*npe*ne*ncu*ncu] x [npe*(npf*nfe)*ne] -> ncu*(npf*nfe)*ncu*(npf*nfe)*ne
void schurMatrixGMinvE(dstype* H, const dstype *G, const dstype *MinvE, dstype scalar, const int npe, const int ncu12, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int M = ndf*npe;    
    int L = ndf*npe*ne; 
    int P = L*ncu12;   
    int N = ncu12*(npf*nfe)*ncu*(npf*nfe)*ne;    
    Kokkos::parallel_for("schurMatrixGMinvE", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu12; // [0, ncu]
        int q = idx/ncu12; // [0, (npf*nfe)*ncu*(npf*nfe)*ne]                  
        int i = q%ndf;   // [0, npf*nfe]
        int p = q/ndf;   // [0, ncu*(npf*nfe)*ne)        
        int n = p%ncu;   // [0, ncu]
        int s = p/ncu;   // [0, (npf*nfe)*ne]                        
        int j = s%ndf;   // [0, (npf*nfe)]
        int e = s/ndf;   // [0, ne]
        dstype sum = 0;
        for (int k=0; k<npe; k++)
            sum += scalar*G[i + ndf*k + M*e + L*m + P*n]*MinvE[k + npe*j + M*e];
        H[idx] -= sum;        
    });                        
}

void schurMatrixGintMinvE(dstype* H, const dstype *G, const dstype *MinvE, dstype scalar, const int npe, const int ncu12, const int ncu, const int npf, const int nfe, const int ne)
{        
    int ndf = npf*nfe;
    int Q = npe*ndf;
    int M = npf*npe;    
    int L = npf*npe*ne; 
    int P = L*ncu12;   
    int N = ncu12*npf*ncu*ndf*ne;    
    Kokkos::parallel_for("schurMatrixGMinvE", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu12; // [0, ncu12]
        int q = idx/ncu12; // [0, npf*ncu*npf*nfe*ne]                  
        int i = q%npf;   // [0, npf]
        int p = q/npf;   // [0, ncu*npf*nfe*ne)        
        int n = p%ncu;   // [0, ncu]
        int s = p/ncu;   // [0, npf*nfe*ne]                        
        int j = s%ndf;   // [0, npf*nfe]
        int e = s/ndf;   // [0, ne]
        dstype sum = 0;
        for (int k=0; k<npe; k++)
            sum += scalar*G[i + npf*k + M*e + L*m + P*n]*MinvE[k + npe*j + Q*e];
        H[idx] -= sum;        
    });                        
}

// npe*ne*ncu-> npe*ncu*ne
void schurVectorRu(dstype* Ru, const dstype* Rutmp,  const int npe, const int ncu, const int ne)
{            
    int L = npe*ne;    
    int N = npe*ncu*ne;    
    Kokkos::parallel_for("schurVectorRu", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%npe; // [0, npe)
        int q = idx/npe; // [0, ncu*ne)                  
        int m = q%ncu;   // [0, ncu]
        int e = q/ncu;   // [0, ne)        
        Ru[idx] = Rutmp[i +  npe*e + L*m];
    });                        
}

// ndf*ne*ncu-> ncu*ndf*ne
void schurVectorRh(dstype* Rh, const dstype* Rhtmp,  const int ndf, const int ncu, const int ne)
{            
    int L = ndf*ne;    
    int N = ndf*ncu*ne;    
    Kokkos::parallel_for("schurVectorRh", N, KOKKOS_LAMBDA(const size_t idx) {
        int m = idx%ncu; // [0, ncu)
        int q = idx/ncu; // [0, ndf*ne)                  
        int i = q%ndf;   // [0, ndf]
        int e = q/ndf;   // [0, ne)        
        Rh[idx] = Rhtmp[i +  ndf*e + L*m];
    });                        
}

void PutFaceNodes(dstype* udg, const dstype* uh, const int* rowe2f1, const int* cole2f1, const int* ent2ind1,
        const int* rowe2f2, const int* cole2f2, const int* ent2ind2, const int npf, const int npe, const int nc, const int e1, const int e2, const int opts)
{
    int ne = e2-e1;
    int K = npf*nc;
    int M = npe*nc;
    int N = M*ne;        
    int I = M*e1;
    if (opts==0) {
        Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            int q, p, s;
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];
            for (j=0; j<n; j++) {
                q = cole2f1[rowe2f1[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] -=  uh[p+npf*m+K*s];                 
            }            
            
            i = ent2ind2[l+npe*k];
            e = (i > 0) ? i : 0;
            n = rowe2f2[i+1] - rowe2f2[e];
            for (j=0; j<n; j++) {
                q = cole2f2[rowe2f2[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]
                udg[I+idx] +=  uh[p+npf*m+K*s];
            }
        });
    }
    else {
        Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];
            for (j=0; j<n; j++) {
                int q = cole2f1[rowe2f1[i]+j];
                int p = q%npf;          // [1, npf]
                int s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] -=  uh[p+npf*m+K*s];
            }            
        });
    }
}

void ApplyXx4(dstype* rg, const dstype* sg, const dstype* fg, const dstype* Xx, const dstype* jac, const int nge, const int nd, const int ncu, const int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*(nd+1);
    int J = I*ncu;
    Kokkos::parallel_for("ApplyXx4", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%M;   // [1, nge*ne]         
        int k = idx/M;   // [1, ncu]
        int g = i%nge;   // [1, nge]
        int e = i/nge;   // [1, ne]
        int ge = g+nge*e;
        int ke = I*k+J*e;
        rg[g+nge*0+ke] = sg[idx]*jac[i];   
        for (int m=0; m<nd; m++) {
            int gem = ge + P*m;
            dstype sum = fg[idx+N*0]*Xx[gem+M*0];
            for (int j=1; j<nd; j++)
                sum += fg[idx+N*j]*Xx[gem+M*j];
            rg[g+nge*(m+1)+ke] = sum;  // nge *(nd+1) * ncu * ne
        }        
    });    
}

void ApplyXxJac(dstype* rg, const dstype* sg, const dstype* fg, const dstype* Xx, const dstype* jac, const int nge, const int nd, const int ncu, const int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*(nd+1);
    int J = I*ne;
    Kokkos::parallel_for("ApplyXx4", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%M;   // [1, nge*ne]         
        int k = idx/M;   // [1, ncu]
        int g = i%nge;   // [1, nge]
        int e = i/nge;   // [1, ne]
        int ge = g+nge*e;
        int ke = I*e + J*k; 
        rg[g+nge*0+ke] = sg[idx]*jac[i];   
        for (int m=0; m<nd; m++) {
            int gem = ge + P*m;
            dstype sum = fg[idx+N*0]*Xx[gem+M*0];
            for (int j=1; j<nd; j++)
                sum += fg[idx+N*j]*Xx[gem+M*j];
            rg[g+nge*(m+1)+ke] = sum;  // nge *(nd+1) * ncu * ne
        }        
    });    
}

void ApplyXxJac(dstype* rg, const dstype* sg_udg, const dstype* fg_udg, const dstype* Xx, const dstype* jac, const int nge, const int nd, const int ncu, const int nc, const int ne)
{
    int M = nge*ne;
    int L = M*ncu;
    int N = M*ncu*nc;
    int P = M*nd;
    int I = nge*(nd+1);
    int J = I*ne;
    int K = J*ncu;
    Kokkos::parallel_for("ApplyXx4", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%M;   // [1, nge*ne]         
        int t = idx/M;   // [1, ncu*nc]
        int k = t%ncu;   // [1, ncu]
        int n = t/ncu;   // [1, nc]
        int g = i%nge;   // [1, nge]
        int e = i/nge;   // [1, ne]
        int ge = g+nge*e;
        int kne = I*e + J*k + K*n; // nge*(nd+1)*e + nge*(nd+1)*ne*k + nge*(nd+1)*ne*nc*n
        rg[g+nge*0+kne] = sg_udg[idx]*jac[i];   
        for (int m=0; m<nd; m++) {
            int gem = ge + P*m;
            int gekn = ge + M*k + L*nd*n;
            dstype sum = fg_udg[gekn + L*0]*Xx[gem+M*0];
            for (int j=1; j<nd; j++)
                sum += fg_udg[gekn + L*j]*Xx[gem+M*j];
            rg[g+nge*(m+1)+kne] = sum;  // nge *(nd+1) * ncu * nc * ne
        }        
    });    
}

// (nge * ne * ncu) x (nge * ne) x (npe * nge) -> (npe * ncu * ne) 
void RuSource(dstype* Ru, const dstype* sg, const dstype* jac, const dstype* testshap, const int nge, const int npe, const int ncu, const int ne)
{
    int N = npe*ncu*ne;
    int M = nge*ne;
    Kokkos::parallel_for("RuSource", N, KOKKOS_LAMBDA(const size_t idx) {
        int p = idx%npe; // [1, npe]         
        int k = idx/npe;   
        int n = k%ncu;   // [1, ncu]
        int e = k/ncu;   // [1, ne]
        dstype sum = 0.0;        
        for (int j=0; j<nge; j++) {            
            sum += sg[j + nge*e + M*n]*jac[j+nge*e]*testshap[p+npe*j];
        }
        Ru[idx] = sum;        
    });    
}

// (nge * ne * ncu * nd) x (nge * ne * nd * nd) x (npe * nge * nd) -> (npe * ncu * ne) 
void RuFlux(dstype* Ru, const dstype* fg, const dstype* Xx, const dstype* testshapderiv, const int nge, const int npe, const int ncu,  int nd, const int ne)
{
    int N = npe*ncu*ne;
    int M = nge*ne;
    int K = nge*ne*ncu;
    int P = nge*ne*nd;
    int Q = npe*nge;
    Kokkos::parallel_for("RuFlux", N, KOKKOS_LAMBDA(const size_t idx) {
        int p = idx%npe; // [1, npe]         
        int k = idx/npe;   
        int n = k%ncu;   // [1, ncu]
        int e = k/ncu;   // [1, ne]
        dstype sum = 0.0;        
        for (int j=0; j<nge; j++) {            
          for (int m=0; m<nd; m++) {
            for (int l=0; l<nd; l++) {
              sum += fg[j + nge*e + M*n + K*l]*Xx[j + nge*e + M*l + P*m]*testshapderiv[p + npe*j + Q*m];
            }
          }
        }
        Ru[idx] += sum;        
    });    
}

// (nge * ne * ncu * nc) x (nge * ne) x (npe * nge) x (nge * npe) -> (npe * ncu * npe * nc * ne) 
void JuSource(dstype* Ju, const dstype* sg_udg, const dstype* jac, const dstype* testshap, const dstype* trialshap, const int nge, const int npe, const int ncu, const int nc, const int ne)
{
    int N = npe*ncu*ne;
    int M = nge*ne;
    Kokkos::parallel_for("RuSource", N, KOKKOS_LAMBDA(const size_t idx) {
        int p = idx%npe; // [1, npe]         
        int tm = idx/npe;  
        int n = tm%ncu;   // [1, ncu]
        int tn = tm/ncu;   //   
        int q = tn%npe;   // [1, npe]
        int tk = tn/npe;   //
        int m = tk%nc;   // [1, nc]
        int e = tk/nc;   // [1, ne]
        dstype sum = 0.0;        
        for (int j=0; j<nge; j++) {            
            sum += sg_udg[j + nge*e + M*n + M*ncu*m]*jac[j+nge*e]*testshap[p+npe*j]*trialshap[j+nge*q];
        }
        Ju[idx] = sum;        
    });    
}

void ApplyXx3(dstype* sg, const dstype* ug, const dstype* Xx, const int nge, const int nd, const int ncu, const int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*nd;
    int J = I*ncu;
    int K = J*nd;
    Kokkos::parallel_for("ApplyXx3", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%M;       // [1, nge*ne]         
        int k = idx/M;   // [1, ncu]
        int g = i%nge;       // [1, nge]
        int e = i/nge;   // [1, ne]
        for (int m=0; m<nd; m++)
            for (int j=0; j<nd; j++)
                sg[g+nge*j+I*k+J*m+K*e] = ug[idx]*Xx[g+nge*e+M*m+P*j]; // nge*nd*ncu*nd*ne
    });    
}

void ApplyJacNormal(dstype* fqg, const dstype* uhg, const dstype* nlg, const dstype* jac, const int nga, const int ncu, const int nd, const int ngf)
{
    int N = nga*ncu;
    int M = ngf*ncu;
    int P = M*nd; 
    Kokkos::parallel_for("ApplyJacNormal", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%nga;   
        int m = idx/nga; // [1, ncu]
        int g = i%ngf;       // [1, ngf]
        int e = i/ngf;   // [1, nf]
        for (int j=0; j<nd; j++)
            fqg[g+ngf*m+M*j+P*e] = uhg[idx]*nlg[i+nga*j]*jac[i];
    });            
}

void ApplyJacFhat(dstype* sg, const dstype* fhg, const dstype* jac, const int nga, const int ncu, const int ngf)
{
    int M = ngf*ncu;
    int N = nga*ncu;
    Kokkos::parallel_for("ApplyJac", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%nga;   
        int m = idx/nga; // [1, ncu]
        int g = i%ngf;       // [1, ngf]
        int e = i/ngf;   // [1, nf]                
        sg[g+ngf*m+M*e] = fhg[idx]*jac[i];
    });
}

void ApplyDtcoef(dstype* sg_u, dstype* fg, dstype dtcoeff, const int nga, const int ncu)
{    
    int N = nga*ncu;
    Kokkos::parallel_for("ApplyJac", N, KOKKOS_LAMBDA(const size_t idx) {
        int i = idx%nga; // [1, nga]  
        int m = idx/nga; // [1, ncu]        
        sg_u[i + nga*m + nga*ncu*m] += dtcoeff*fg[i + nga*m];
    });
}

void ApplyJac(dstype* R, const dstype* jac, const int M, const int N)
{
    // M = npe*ncr
    // N = npe*ncr*ne
    // R:   [npe*ncr*ne]
    // jac: [ne]    
    Kokkos::parallel_for("ApplyJac", N, KOKKOS_LAMBDA(const size_t n) {            
        int m = n%M;       // [1, npe*ncr] 
        int i = (n-m)/M;   // [1, ne] 
        R[n] = R[n]*(jac[i]);
    });                        
}

void ApplyJacInv(dstype* R, const dstype* jac, const int M, const int N)
{
    // M = npe*ncr
    // N = npe*ncr*ne
    // R:   [npe*ncr*ne]
    // jac: [ne]    
    Kokkos::parallel_for("ApplyJacInv", N, KOKKOS_LAMBDA(const size_t n) {            
        int m = n%M;       // [1, npe*ncr] 
        int i = (n-m)/M;   // [1, ne] 
        R[n] = R[n]/(jac[i]);
    });                        
}

void ShapJac(dstype* shapjac, const dstype* shapegt, const dstype* jac, const int nge, const int npe, const int ne)
{
    int M = nge*npe;
    int N = nge*npe*ne;
    Kokkos::parallel_for("ShapJac", N, KOKKOS_LAMBDA(const size_t i) {
        int l = i%M;       // [1, nge*npe]           
        int k = (i-l)/M;   // [1, ne] 
        int n = l%nge;     // [1,nge]
        int m = (l-n)/nge; // [1, npe]
        shapjac[i] = shapegt[n+nge*m]*jac[n+nge*k];
    });    
}

void ElemGeom1D(dstype* jac, dstype* Xx, const dstype* Jg, const int ngv)
{
    Kokkos::parallel_for("ElemGeom1D", ngv, KOKKOS_LAMBDA(const size_t i) {
        jac[i] = Jg[i];
        Xx[i] = 1.0;
    });
}

void ElemGeom2D(dstype* jac, dstype* Xx11, dstype* Xx12, dstype* Xx21, dstype* Xx22,
                 const dstype* Jg11, const dstype* Jg12, const dstype* Jg21, const dstype* Jg22, const int ngv)
{
    Kokkos::parallel_for("ElemGeom2D", ngv, KOKKOS_LAMBDA(const size_t i) {
        jac[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
        Xx11[i] = Jg22[i]; // dxi/dx
        Xx21[i] = -Jg21[i]; //dxi/dy 
        Xx12[i] = -Jg12[i]; //deta/dx
        Xx22[i] = Jg11[i]; //deta/dy
    });
}

void ElemGeom3D(dstype* jac, dstype* Xx11, dstype* Xx12, dstype* Xx13, dstype* Xx21, 
                dstype* Xx22, dstype* Xx23, dstype* Xx31, dstype* Xx32, dstype* Xx33,
                const dstype* Jg11, const dstype* Jg12, const dstype* Jg13, const dstype* Jg21, const dstype* Jg22, 
                const dstype* Jg23, const dstype* Jg31, const dstype* Jg32, const dstype* Jg33, const int ngv)
{
    Kokkos::parallel_for("ElemGeom3D", ngv, KOKKOS_LAMBDA(const size_t i) {
        jac[i] = Jg11[i]*Jg22[i]*Jg33[i] - Jg11[i]*Jg32[i]*Jg23[i] +
                 Jg21[i]*Jg32[i]*Jg13[i] - Jg21[i]*Jg12[i]*Jg33[i] +
                 Jg31[i]*Jg12[i]*Jg23[i] - Jg31[i]*Jg22[i]*Jg13[i];
        Xx11[i] = Jg22[i]*Jg33[i] - Jg23[i]*Jg32[i];
        Xx21[i] = Jg23[i]*Jg31[i] - Jg21[i]*Jg33[i];
        Xx31[i] = Jg21[i]*Jg32[i] - Jg22[i]*Jg31[i];
        Xx12[i] = Jg13[i]*Jg32[i] - Jg12[i]*Jg33[i];
        Xx22[i] = Jg11[i]*Jg33[i] - Jg13[i]*Jg31[i];
        Xx32[i] = Jg12[i]*Jg31[i] - Jg11[i]*Jg32[i];
        Xx13[i] = Jg12[i]*Jg23[i] - Jg13[i]*Jg22[i];
        Xx23[i] = Jg13[i]*Jg21[i] - Jg11[i]*Jg23[i];
        Xx33[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];        
    });
}

void FaceGeom1D(dstype* jacg, dstype* nlg, const dstype* Jg, const int na)
{       
    Kokkos::parallel_for("FaceGeom1D", na, KOKKOS_LAMBDA(const size_t i) {
        jacg[i] = 1.0;
        nlg[i] = 1.0;
    });
}

void FixNormal1D(dstype* nlg, const int* facecon, const int na)
{       
    Kokkos::parallel_for("FixNormal1D", na, KOKKOS_LAMBDA(const size_t i) {
        if (facecon[2*i]==0)
            nlg[i] = -1.0;
    });
}

void FixNormal1D(dstype* nlg, const int na)
{       
    Kokkos::parallel_for("FixNormal1D", na, KOKKOS_LAMBDA(const size_t i) {
        if (i%2 == 0)
            nlg[i] = -1.0;
    });
}

void FaceGeom2D(dstype* jacg, dstype* nlg, const dstype* Jg, const int na)
{       
    Kokkos::parallel_for("FaceGeom2D", na, KOKKOS_LAMBDA(const size_t i) {
        int j = i+na;
        jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
        nlg[i] = Jg[j]/jacg[i];
        nlg[j] = -Jg[i]/jacg[i];
    });
}

void FaceGeom3D(dstype* jacg, dstype* nlg, const dstype* Jg, const int na)
{       
    int n11 = 0;
    int n21 = na;
    int n31 = 2*na;
    int n12 = 3*na;
    int n22 = 4*na;
    int n32 = 5*na;
    Kokkos::parallel_for("FaceGeom3D", na, KOKKOS_LAMBDA(const size_t i) {
        int j = i+na;
        int k = i+2*na;
        nlg[i] = Jg[i+n21]*Jg[i+n32] - Jg[i+n31]*Jg[i+n22];
        nlg[j] = Jg[i+n31]*Jg[i+n12] - Jg[i+n11]*Jg[i+n32];
        nlg[k] = Jg[i+n11]*Jg[i+n22] - Jg[i+n21]*Jg[i+n12];
        jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
        nlg[i] = nlg[i]/jacg[i];
        nlg[j] = nlg[j]/jacg[i];
        nlg[k] = nlg[k]/jacg[i];
    });
}

#endif  


// *************************************************************************************** //

// void AverageFlux(view_1d fg, const int N)
// {	
//     Kokkos::parallel_for("AverageFlux", N, KOKKOS_LAMBDA(const size_t tid) {
//         fg(tid+N) = 0.5*(fg(tid) + fg(tid+N));            
//     });
// }

// void AverageFluxDotNormal(view_1d fg, const view_1d nl, const int N, const int M, const int numPoints, const int nd)
// {	
//     Kokkos::parallel_for("AverageFluxDotNormal", M, KOKKOS_LAMBDA(const size_t tid) {
//         int i = tid%numPoints;                
//         dstype sum = fg(N + 0*M + tid) * nl(i + 0 * numPoints);   
//         for (int l = 1; l < nd; l++)
//             sum += fg(N + l*M + tid) * nl(i + l * numPoints);
//         fg(tid) = sum;     
//     });
// }

// void AddStabilization1(view_1d fg, const view_1d ug1, const view_1d ug2, const view_1d tau, const int M)
// {	
//     Kokkos::parallel_for("AddStabilization1", M, KOKKOS_LAMBDA(const size_t tid) {
//         fg(tid) += tau(0) * (ug1(tid) - ug2(tid));        
//     });
// }

// void AddStabilization2(view_1d fg, const view_1d ug1, const view_1d ug2, const view_1d tau, const int M, const int numPoints)
// {	
//     Kokkos::parallel_for("AddStabilization2", M, KOKKOS_LAMBDA(const size_t tid) {
//         int i = tid%numPoints;   
//         int j = (tid-i)/numPoints;
//         fg(tid) += tau(j) * (ug1(tid) - ug2(tid));     
//     });
// }

// void AddStabilization3(view_1d fg, const view_1d ug1, const view_1d ug2, const view_1d tau, const int M, const int numPoints, const int ncu)
// {	
//     Kokkos::parallel_for("AddStabilization3", M, KOKKOS_LAMBDA(const size_t tid) {
//         int i = tid%numPoints;   
//         int j = (tid-i)/numPoints;
//         for (int k=0; k<ncu; k++) {
//             int nm = k * ncu + j;
//             int nk = k * numPoints + i;
//             fg(tid) += tau(nm) * (ug1(nk) - ug2(nk));
//         }
//     });
// }

// void GetArrayAtIndex(view_1d y, const view_1d x, const view_1i ind, const int n)
// {    
//     Kokkos::parallel_for("GetArrayAtIndex", n, KOKKOS_LAMBDA(const size_t i) {
//         y(i) = x(ind(i));
//     });
// }

// void PutArrayAtIndex(view_1d y, const view_1d x, const view_1i ind, const int n)
// {    
//     Kokkos::parallel_for("PutArrayAtIndex", n, KOKKOS_LAMBDA(const size_t i) {
//         y(ind(i)) = x(i);
//     });
// }

// void ArrayCopy(view_1d y, const view_1d x, const int n)
// {    
//     Kokkos::parallel_for("GetArrayAtIndex", n, KOKKOS_LAMBDA(const size_t i) {
//         y(i) = x(i);
//     });
// }

// void ArraySetValue(view_1d y, const dstype a, const int n)
// {    
//     Kokkos::parallel_for("ArraySetValue", n, KOKKOS_LAMBDA(const size_t i) {
//         y(i) = a;
//     });
// }

// void ArrayMultiplyScalar(view_1d y, const dstype a, const int n)
// {    
//     Kokkos::parallel_for("ArrayMultiplyScalar", n, KOKKOS_LAMBDA(const size_t i) {
//         y(i) = a*y(i);
//     });
// }

// void ArrayAXPB(view_1d y, view_1d x, const dstype a, const dstype b, const int N) 
// {
//     Kokkos::parallel_for("ArrayAXPBY", N, KOKKOS_LAMBDA(const size_t idx) {
//         y(idx) = a * x(idx) + b;
//     });
// }

// void ArrayAXPBY(view_1d y, view_1d x, const view_1d z, const dstype a, const dstype b, const int N) 
// {
//     Kokkos::parallel_for("ArrayAXPBY", N, KOKKOS_LAMBDA(const size_t idx) {
//         y(idx) = a * x(idx) + b * z(idx);
//     });
// }

// void ArrayAXY(view_1d y, view_1d x, const view_1d z, const dstype a, const int N) 
// {
//     Kokkos::parallel_for("ArrayAXY", N, KOKKOS_LAMBDA(const size_t idx) {
//         y(idx) = a * x(idx) * z(idx);
//     });
// }

// void ArrayExtract(view_1d un, const view_1d u, const int I, const int J, const int K, 
//         const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
// {        
//     int ni = i2-i1;
//     int nj = j2-j1;
//     int nk = k2-k1;    
//     int M = ni*nj;
//     int N = M*nk;
//     int Q = I*J;
//     Kokkos::parallel_for("ArrayExtract", N, KOKKOS_LAMBDA(const size_t idx) {
//         int l = idx%M;        
//         int i = l%ni+i1;
//         int j = (l-i)/ni+j1;
//         int k = (idx-l)/M+k1;
//         un(idx) = u(i+I*j+Q*k);
//     });
// }

// void ArrayInsert(view_1d u, const view_1d un, const int I, const int J, const int K, 
//         const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
// {        
//     int ni = i2-i1;
//     int nj = j2-j1;
//     int nk = k2-k1;    
//     int M = ni*nj;
//     int N = M*nk;
//     int Q = I*J;
//     Kokkos::parallel_for("ArrayInsert", N, KOKKOS_LAMBDA(const size_t idx) {
//         int l = idx%M;        
//         int i = l%ni+i1;
//         int j = (l-i)/ni+j1;
//         int k = (idx-l)/M+k1;
//         u(i+I*j+Q*k) = un(idx);        
//     });
// }

// void ArrayGemmBatch(view_1d C, const view_1d A, const view_1d B, const int I, const int J, const int K, const int S)
// {        
//     // C[I*J*S] = A[I*K*S] x B[K*J*S]
//     int M = I*J;
//     int N = M*S;
//     Kokkos::parallel_for("ArrayGemmBatch", N, KOKKOS_LAMBDA(const size_t idx) {
//         int l = idx%M;        
//         int i = l%I;
//         int j = (l-i)/I;
//         int s = (idx-l)/M;
//         dstype sum = 0.0;
//         for (int k=0; k<K; k++)
//             sum += A(i+I*k+I*K*s)*B(k+K*j+K*J*s);
//         C(idx) = sum;    
//     });
// }

// void ArrayGemmBatch1(view_1d C, const view_1d A, const view_1d B, const int I, const int J, const int K, const int S)
// {        
//     // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
//     int M = I*J;
//     int N = M*S;
//     int Q = I*K;
//     int P = K*J;
//     Kokkos::parallel_for("ArrayGemmBatch1", N, KOKKOS_LAMBDA(const size_t idx) {
//         int l = idx%M;        
//         int i = l%I;
//         int j = (l-i)/I;
//         int s = (idx-l)/M;
//         dstype sum = C(idx);
//         for (int k=0; k<K; k++)
//             sum += A(i+I*k+Q*s)*B(k+K*j+P*s);
//        C(idx) = sum;     
//     });
// }

// void ArrayDG2CG(view_1d ucg, const view_1d udg, const view_1i cgent2dgent, const view_1i rowent2elem, const int nent)
// {        
//     Kokkos::parallel_for("ArrayDG2CG", nent, KOKKOS_LAMBDA(const size_t i) {
//         dstype sum = 0.0;
//         int nelem = rowent2elem(i+1)-rowent2elem(i);
//         dstype fac = 1.0/((dstype) nelem);
//         for (int k=0; k<nelem; k++)
//             sum += udg(cgent2dgent(rowent2elem(i)+k)); 
//         ucg(i) = sum*fac;
//     });
// }

// void ArrayDG2CG2(view_1d ucg, const view_1d udg, const view_1i colent2elem, const view_1i rowent2elem, const int nent, const int npe)
// {        
//     Kokkos::parallel_for("ArrayDG2CG2", nent, KOKKOS_LAMBDA(const size_t i) {        
//         int nelem = rowent2elem(i+1)-rowent2elem(i);
//         dstype fac = 1.0/((dstype) (nelem*npe));
//         dstype sum = 0.0;
//         for (int k=0; k<nelem; k++) {
//             int e = colent2elem(rowent2elem(i)+k);
//             for (int j=0; j<npe; j++)
//                 sum += udg(j+npe*e); 
//         }
//         ucg(i) = sum*fac;
//     });
// }

// void ArrayMatrixMultiplication(view_1d C, const view_1d A, const view_1d B, const int I, const int J, const int K)
// {        
//     // C[I*J] = A[I*K] x B[K*J]
//     int N = I*J;    
//     Kokkos::parallel_for("ArrayMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%I;   //   [1, I]
//         int j = idx/I; // [1, J]        
//         dstype sum = 0.0;
//         for (int k=0; k<K; k++)
//             sum += A(i + I*k)*B(k + K*j);
//         C(idx) = sum;    
//     });
// }

// void ArrayMatrixMultiplication1(view_1d C, const view_1d A, const view_1d B, const int I, const int J, const int K)
// {        
//     // C[I*J] += A[I*K] x B[K*J]
//     int N = I*J;    
//     Kokkos::parallel_for("ArrayMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%I;   //   [1, I]
//         int j = idx/I; // [1, J]        
//         dstype sum = 0.0;
//         for (int k=0; k<K; k++)
//             sum += A(i + I*k)*B(k + K*j);
//         C(idx) += sum;    
//     });
// }

// void ArrayEosInverseMatrix11(view_1d A, const int npe, const int ncw, const int ne)
// {            
//     int N = npe*ne;
//     Kokkos::parallel_for("ArrayEosInverseMatrix11", N, KOKKOS_LAMBDA(const size_t i) {
//         A(i) = 1.0/A(i);    
//     });
// }

// void ArrayEosInverseMatrix22(view_1d A, const int npe, const int ncw, const int ne)
// {        
//     int N = npe*ne;
//     int M = npe*ncw*ncw;
//     Kokkos::parallel_for("ArrayEosInverseMatrix22", N, KOKKOS_LAMBDA(const size_t i) {
//         int j = i%npe;    // [1, npe]
//         int k = (i-j)/npe; //[1, ne]
//         int jk = j + M*k;
//         dstype a11 = A(jk + npe*0);
//         dstype a21 = A(jk + npe*1);
//         dstype a12 = A(jk + npe*2);
//         dstype a22 = A(jk + npe*3);
//         dstype detA = (a11*a22- a12*a21);      
//         A(jk + npe*0) = a22/detA;
//         A(jk + npe*1) = -a21/detA;
//         A(jk + npe*2) = -a12/detA;
//         A(jk + npe*3) = a11/detA;
//     });
// }

// void ArrayEosInverseMatrix33(view_1d A, const int npe, const int ncw, const int ne)
// {            
//     int N = npe*ne;
//     int M = npe*ncw*ncw;
//     Kokkos::parallel_for("ArrayEosInverseMatrix33", N, KOKKOS_LAMBDA(const size_t i) {
//         int j = i%npe;    // [1, npe]
//         int k = (i-j)/npe; //[1, ne]
//         int jk = j + M*k;
//         dstype a11 = A(jk + npe*0);
//         dstype a21 = A(jk + npe*1);
//         dstype a31 = A(jk + npe*2);
//         dstype a12 = A(jk + npe*3);
//         dstype a22 = A(jk + npe*4);
//         dstype a32 = A(jk + npe*5);
//         dstype a13 = A(jk + npe*6);
//         dstype a23 = A(jk + npe*7);
//         dstype a33 = A(jk + npe*8);        
//         dstype detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
      
//         A(jk + npe*0) = (a22*a33 - a23*a32)/detA;
//         A(jk + npe*1) = (a23*a31 - a21*a33)/detA;
//         A(jk + npe*2) = (a21*a32 - a22*a31)/detA;
//         A(jk + npe*3) = (a13*a32 - a12*a33)/detA;
//         A(jk + npe*4) = (a11*a33 - a13*a31)/detA;
//         A(jk + npe*5) = (a12*a31 - a11*a32)/detA;
//         A(jk + npe*6) = (a12*a23 - a13*a22)/detA;
//         A(jk + npe*7) = (a13*a21 - a11*a23)/detA;
//         A(jk + npe*8) = (a11*a22 - a12*a21)/detA;
//     });
// }

// void ArrayEosMatrixMultiplication(view_1d C, const view_1d A, const view_1d B, const int npe, const int ncw, const int ne, const int ncu)
// {        
//     // C[npe*ncw*ncu*ne] = A[npe*ncw*ncw*ne] x B[npe*ncw*ncu*ne]
//     int N = npe*ne;
//     int K = npe*ncw;
//     int P = K*ncu;
//     int Q = K*ncw;    
//     Kokkos::parallel_for("ArrayEosMatrixMultiplication", N, KOKKOS_LAMBDA(const size_t i) {
//         int j = i%npe;    // [1, npe]
//         int k = (i-j)/npe; //[1, ne]        
//         for (int b=0; b<ncu; b++)
//           for (int a=0; a<ncw; a++) {
//             dstype sum = 0.0;
//             for (int m=0; m<ncw; m++)
//               sum += A(j + npe*a + K*m + Q*k)*B(j + npe*m + K*b + P*k);
//             C(j + npe*a + K*b + P*k) = sum;  
//           }        
//     });
// }

// void GetElemNodes(view_1d unView, const view_1d uView, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
// {
//     int nn = np * (e2 - e1);
//     int ncu = nc2 - nc1;
//     int N = nn * ncu;
//     int K = np * nc;

//     Kokkos::parallel_for("GetElemNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx % nn;  // [0, np*(e2-e1)]
//         int j = idx / nn;  // [0, ncu]
//         int k = i % np;    // [0, np]
//         int e = i / np + e1;
//         unView(idx) = uView(k + (j + nc1) * np + e * K);
//     });
// }

// void PutElemNodes(view_1d u, const view_1d un, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
// {
//     int nn = np * (e2 - e1);
//     int ncu = nc2 - nc1;
//     int N = nn * ncu;
//     int K = np * nc;
//     Kokkos::parallel_for("PutElemNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx % nn;  // [0, np*(e2-e1)]
//         int j = idx / nn;  // [0, ncu]
//         int k = i % np;    // [0, np]
//         int e = i / np + e1;
//         u(k + (j + nc1) * np + e * K) = un(idx);
//     });
// }

// void GetFaceNodes(view_1d uh, const view_1d udg, const view_1i facecon, const int npf, const int ncu, const int npe, const int nc, const int f1, const int f2, const int opts) 
//  {
//     int nf = f2-f1;
//     int ndf = npf*nf;
//     int N = ndf*ncu;
//     int M = npe*nc;
//     Kokkos::parallel_for("GetFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%ndf; 
//         int j = idx/ndf; // [0, ncu)
//         int m = npf*f1+i;
//         int k1 = facecon(2*m);
//         int k2 = facecon(2*m+1);
//         int m1 = k1%npe; // [0, npe)
//         int m2 = k2%npe; // [0, npe)
//         int n1 = (k1-m1)/npe; // [0, ne)
//         int n2 = (k2-m2)/npe; // [0, ne)              
//         // uh npf*nf*ncu 
//         if (opts==0) 
//             uh(idx) = 0.5*(udg(m1+j*npe+n1*M)+udg(m2+j*npe+n2*M));
//         else if (opts==1) 
//             uh(idx) = udg(m1+j*npe+n1*M);
//         else if (opts==2) 
//             uh(idx) = udg(m2+j*npe+n2*M);
//     });
//  }

// void PutFaceNodes(view_1d udg, const view_1d uh, const view_1i facecon, const int npf, const int ncu, const int npe, const int nc, const int f1, const int f2)
// {
//     int nf = f2-f1;
//     int ndf = npf*nf;
//     int N = ndf*ncu;
//     int M = npe*nc;
//     Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//         int p = idx%npf; // [0, npf)
//         int q = idx/npf; // [0, ncu*nf)        
//         int j = q%ncu;  // [0, ncu]
//         int f = q/ncu; // [0, nf)
//         int i = p + npf*f;
//         int m = npf*f1+i;
//         int k1 = facecon(2*m);
//         int k2 = facecon(2*m+1);
//         int m1 = k1%npe; // [0, npe)
//         int m2 = k2%npe; // [0, npe)
//         int n1 = (k1-m1)/npe; // [0, ne)
//         int n2 = (k2-m2)/npe; // [0, ne)                          
//         dstype uhx = uh(idx); // npf*ncu*nf               
//         if (k1==k2) { // boundary face 
//             Kokkos::atomic_sub(&udg(m1+j*npe+n1*M), uhx);            
//         }
//         else { // interior face
//             Kokkos::atomic_sub(&udg(m1+j*npe+n1*M), uhx);
//             Kokkos::atomic_add(&udg(m2+j*npe+n2*M), uhx);
//         }        
//     });                        
// }

// void PutFaceNodes(view_1d udg, const view_1d uh, const view_1i rowe2f1, const view_1i cole2f1, const view_1i ent2ind1,
//         const view_1i rowe2f2, const view_1i cole2f2, const view_1i ent2ind2, const int npf, const int npe, const int nc, const int e1, const int e2, const int opts)
// {
//     int ne = e2-e1;
//     int K = npf*nc;
//     int M = npe*nc;
//     int N = M*ne;        
//     int I = M*e1;
//     if (opts==0) {
//         Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//             int j = idx%M;              //[1, npe*nc]
//             int k = (idx-j)/M+e1;       //[1, ne]      
//             int l = j%npe;              //[1, npe]
//             int m = (j-l)/npe;          //[1, nc] 
//             int q, p, s;
            
//             int i = ent2ind1(l+npe*k);
//             int e = (i > 0) ? i : 0;
//             int n = rowe2f1(i+1) - rowe2f1(e);
//             for (j=0; j<n; j++) {
//                 q = cole2f1(rowe2f1(i)+j);
//                 p = q%npf;              // [1, npf]
//                 s = (q-p)/npf;          // [1, nf]    
//                 udg(I+idx) -=  uh(p+npf*m+K*s);                 
//             }            
            
//             i = ent2ind2(l+npe*k);
//             e = (i > 0) ? i : 0;
//             n = rowe2f2(i+1) - rowe2f2(e);
//             for (j=0; j<n; j++) {
//                 q = cole2f2(rowe2f2(i)+j);
//                 p = q%npf;              // [1, npf]
//                 s = (q-p)/npf;          // [1, nf]
//                 udg(I+idx) +=  uh(p+npf*m+K*s);
//             }
//         });
//     }
//     else {
//         Kokkos::parallel_for("PutFaceNodes", N, KOKKOS_LAMBDA(const size_t idx) {
//             int j = idx%M;              //[1, npe*nc]
//             int k = (idx-j)/M+e1;       //[1, ne]      
//             int l = j%npe;              //[1, npe]
//             int m = (j-l)/npe;          //[1, nc] 
            
//             int i = ent2ind1(l+npe*k);
//             int e = (i > 0) ? i : 0;
//             int n = rowe2f1(i+1) - rowe2f1(e);
//             for (j=0; j<n; j++) {
//                 int q = cole2f1(rowe2f1(i)+j);
//                 int p = q%npf;          // [1, npf]
//                 int s = (q-p)/npf;          // [1, nf]    
//                 udg(I+idx) -=  uh(p+npf*m+K*s);
//             }            
//         });
//     }
// }

// void ApplyXx4(view_1d rg, const view_1d sg, const view_1d fg, const view_1d Xx, const view_1d jac, const int nge, const int nd, const int ncu, const int ne)
// {
//     int M = nge*ne;
//     int N = M*ncu;
//     int P = M*nd;
//     int I = nge*(nd+1);
//     int J = I*ncu;
//     Kokkos::parallel_for("ApplyXx4", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%M;   // [1, nge*ne]         
//         int k = idx/M;   // [1, ncu]
//         int g = i%nge;   // [1, nge]
//         int e = i/nge;   // [1, ne]
//         int ge = g+nge*e;
//         int ke = I*k+J*e;
//         rg[g+nge*0+ke] = sg[idx]*jac[i];   
//         for (int m=0; m<nd; m++) {
//             int gem = ge + P*m;
//             dstype sum = fg(idx+N*0)*Xx(gem+M*0);
//             for (int j=1; j<nd; j++)
//                 sum += fg(idx+N*j)*Xx(gem+M*j);
//             rg(g+nge*(m+1)+ke) = sum;  // nge *(nd+1) * ncu * ne
//         }        
//     });    
// }

// void ApplyXx3(view_1d sg, const view_1d ug, const view_1d Xx, const int nge, const int nd, const int ncu, const int ne)
// {
//     int M = nge*ne;
//     int N = M*ncu;
//     int P = M*nd;
//     int I = nge*nd;
//     int J = I*ncu;
//     int K = J*nd;
//     Kokkos::parallel_for("ApplyXx3", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%M;       // [1, nge*ne]         
//         int k = (idx)/M;   // [1, ncu]
//         int g = i%nge;       // [1, nge]
//         int e = (i)/nge;   // [1, ne]
//         for (int m=0; m<nd; m++)
//             for (int j=0; j<nd; j++)
//                 sg(g+nge*j+I*k+J*m+K*e) = ug(idx)*Xx(g+nge*e+M*m+P*j); // nge*nd*ncu*nd*ne
//     });    
// }

// void ApplyJacNormal(view_1d fqg, const view_1d uhg, const view_1d nlg, const view_1d jac, const int nga, const int ncu, const int nd, const int ngf)
// {
//     int N = nga*ncu;
//     int M = ngf*ncu;
//     int P = M*nd; 
//     Kokkos::parallel_for("ApplyJacNormal", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%nga;   
//         int m = (idx)/nga; // [1, ncu]
//         int g = i%ngf;       // [1, ngf]
//         int e = (i)/ngf;   // [1, nf]
//         for (int j=0; j<nd; j++)
//             fqg(g+ngf*m+M*j+P*e) = uhg(idx)*nlg(i+nga*j)*jac(i);
//     });            
// }

// void ApplyJacFhat(view_1d sg, const view_1d fhg, const view_1d jac, const int nga, const int ncu, const int ngf)
// {
//     int M = ngf*ncu;
//     int N = nga*ncu;
//     Kokkos::parallel_for("ApplyJacFhat", N, KOKKOS_LAMBDA(const size_t idx) {
//         int i = idx%nga;   
//         int m = idx/nga; // [1, ncu]
//         int g = i%ngf;       // [1, ngf]
//         int e = i/ngf;   // [1, nf]                
//         sg(g+ngf*m+M*e) = fhg(idx)*jac(i);
//     });
// }


// void ApplyJac(view_1d R, const view_1d jac, const int M, const int N)
// {
//     // M = npe*ncr
//     // N = npe*ncr*ne
//     // R:   [npe*ncr*ne]
//     // jac: [ne]    
//     Kokkos::parallel_for("ApplyJac", N, KOKKOS_LAMBDA(const size_t n) {            
//         int m = n%M;       // [1, npe*ncr] 
//         int i = (n-m)/M;   // [1, ne] 
//         R(n) = R(n)*(jac(i));
//     });                        
// }

// void ApplyJacInv(view_1d R, const view_1d jac, const int M, const int N)
// {
//     // M = npe*ncr
//     // N = npe*ncr*ne
//     // R:   [npe*ncr*ne]
//     // jac: [ne]    
//     Kokkos::parallel_for("ApplyJacInv", N, KOKKOS_LAMBDA(const size_t n) {            
//         int m = n%M;       // [1, npe*ncr] 
//         int i = (n-m)/M;   // [1, ne] 
//         R(n) = R(n)/(jac(i));
//     });                        
// }

// void ShapJac(view_1d shapjac, const view_1d shapegt, const view_1d jac, const int nge, const int npe, const int ne)
// {
//     int M = nge*npe;
//     int N = nge*npe*ne;
//     Kokkos::parallel_for("ShapJac", N, KOKKOS_LAMBDA(const size_t i) {
//         int l = i%M;       // [1, nge*npe]           
//         int k = (i-l)/M;   // [1, ne] 
//         int n = l%nge;     // [1,nge]
//         int m = (l-n)/nge; // [1, npe]
//         shapjac(i) = shapegt(n+nge*m)*jac(n+nge*k);
//     });    
// }

// void ElemGeom1D(view_1d jac, view_1d Xx, const view_1d Jg, const int ngv)
// {
//     Kokkos::parallel_for("ElemGeom1D", ngv, KOKKOS_LAMBDA(const size_t i) {
//         jac(i) = Jg(i);
//         Xx(i) = 1.0;
//     });
// }

// void ElemGeom2D(view_1d jac, view_1d Xx11, view_1d Xx12, view_1d Xx21, view_1d Xx22,
//                  const view_1d Jg11, const view_1d Jg12, const view_1d Jg21, const view_1d Jg22, const int ngv)
// {
//     Kokkos::parallel_for("ElemGeom2D", ngv, KOKKOS_LAMBDA(const size_t i) {
//         jac(i) = Jg11(i)*Jg22(i) - Jg12(i)*Jg21(i);
//         Xx11(i) = Jg22(i); // dxi/dx
//         Xx21(i) = -Jg21(i); //dxi/dy 
//         Xx12(i) = -Jg12(i); //deta/dx
//         Xx22(i) = Jg11(i); //deta/dy
//     });
// }

// void ElemGeom3D(view_1d jac, view_1d Xx11, view_1d Xx12, view_1d Xx13, view_1d Xx21, 
//                 view_1d Xx22, view_1d Xx23, view_1d Xx31, view_1d Xx32, view_1d Xx33,
//                 const view_1d Jg11, const view_1d Jg12, const view_1d Jg13, const view_1d Jg21, const view_1d Jg22, 
//                 const view_1d Jg23, const view_1d Jg31, const view_1d Jg32, const view_1d Jg33, const int ngv)
// {
//     Kokkos::parallel_for("ElemGeom3D", ngv, KOKKOS_LAMBDA(const size_t i) {
//         jac(i) = Jg11(i)*Jg22(i)*Jg33(i) - Jg11(i)*Jg32(i)*Jg23(i) +
//                  Jg21(i)*Jg32(i)*Jg13(i) - Jg21(i)*Jg12(i)*Jg33(i) +
//                  Jg31(i)*Jg12(i)*Jg23(i) - Jg31(i)*Jg22(i)*Jg13(i);
//         Xx11(i) = Jg22(i)*Jg33(i) - Jg23(i)*Jg32(i);
//         Xx21(i) = Jg23(i)*Jg31(i) - Jg21(i)*Jg33(i);
//         Xx31(i) = Jg21(i)*Jg32(i) - Jg22(i)*Jg31(i);
//         Xx12(i) = Jg13(i)*Jg32(i) - Jg12(i)*Jg33(i);
//         Xx22(i) = Jg11(i)*Jg33(i) - Jg13(i)*Jg31(i);
//         Xx32(i) = Jg12(i)*Jg31(i) - Jg11(i)*Jg32(i);
//         Xx13(i) = Jg12(i)*Jg23(i) - Jg13(i)*Jg22(i);
//         Xx23(i) = Jg13(i)*Jg21(i) - Jg11(i)*Jg23(i);
//         Xx33(i) = Jg11(i)*Jg22(i) - Jg12(i)*Jg21(i);        
//     });
// }

// void FaceGeom1D(view_1d jacg, view_1d nlg, const view_1d Jg, const int na)
// {       
//     Kokkos::parallel_for("FaceGeom1D", na, KOKKOS_LAMBDA(const size_t i) {
//         jacg(i) = 1.0;
//         nlg(i) = 1.0;
//     });
// }

// void FixNormal1D(view_1d nlg, const view_1i facecon, const int na)
// {       
//     Kokkos::parallel_for("FixNormal1D", na, KOKKOS_LAMBDA(const size_t i) {
//         if (facecon(2*i)==0)
//             nlg(i) = -1.0;
//     });
// }

// void FaceGeom2D(view_1d jacg, view_1d nlg, const view_1d Jg, const int na)
// {       
//     Kokkos::parallel_for("FaceGeom2D", na, KOKKOS_LAMBDA(const size_t i) {
//         int j = i+na;
//         jacg(i) = sqrt(Jg(i)*Jg(i) + Jg(j)*Jg(j));
//         nlg(i) = Jg(j)/jacg(i);
//         nlg(j) = -Jg(i)/jacg(i);
//     });
// }

// void FaceGeom3D(view_1d jacg, view_1d nlg, const view_1d Jg, const int na)
// {       
//     int n11 = 0;
//     int n21 = na;
//     int n31 = 2*na;
//     int n12 = 3*na;
//     int n22 = 4*na;
//     int n32 = 5*na;
//     Kokkos::parallel_for("FaceGeom3D", na, KOKKOS_LAMBDA(const size_t i) {
//         int j = i+na;
//         int k = i+2*na;
//         nlg(i) = Jg(i+n21)*Jg(i+n32) - Jg(i+n31)*Jg(i+n22);
//         nlg(j) = Jg(i+n31)*Jg(i+n12) - Jg(i+n11)*Jg(i+n32);
//         nlg(k) = Jg(i+n11)*Jg(i+n22) - Jg(i+n21)*Jg(i+n12);
//         jacg(i) = sqrt(nlg(i)*nlg(i) + nlg(j)*nlg(j) + nlg(k)*nlg(k));
//         nlg(i) = nlg(i)/jacg(i);
//         nlg(j) = nlg(j)/jacg(i);
//         nlg(k) = nlg(k)/jacg(i);
//     });
// }
