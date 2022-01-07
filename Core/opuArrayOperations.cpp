#ifndef __OPUARRAYOPERATIONS
#define __OPUARRAYOPERATIONS

template <typename T> void opuGetArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    for (int i = 0; i<n; i++)    
        y[i] = x[ind[i]];
}

template <typename T> void opuPutArrayAtIndex(T *y, T *x, int *ind, int n)
{    
    for (int i = 0; i<n; i++)    
        y[ind[i]] = x[i];
}

template <typename T> void opuArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{    
    for (int i = 0; i<n; i++)    
        y[ind[i]] += x[i];
}

template <typename T> void opuArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{    
    for (int i = 0; i<n; i++)    
        y[ind[i]] -= x[i];
}

template <typename T> void opuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{    
    for (int i = 0; i<n; i++)    
        y[ind[i]] = y[ind[i]] + a*x[i];
}

template <typename T> void opuArrayAverageAtIndex(T *y, T *x, int *ind1, int *ind2, int n)
{    
    for (int i = 0; i<n; i++)    
        y[i] = ((T) 0.5)*(x[ind1[i]] + x[ind2[i]]);
}

template <typename T> T opuArrayMin(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

template <typename T> T opuArrayMax(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

template <typename T> T opuArraySum(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    return b;
}

template <typename T> T opuArraySquareSum(T *a, int n)
{
    T b = a[0]*a[0];
    for (int i=1; i<n; i++)
            b = b + a[i]*a[i];    
    return b;
}

template <typename T> T opuArrayMean(T *a, int n)
{
    T b = a[0];
    for (int i=1; i<n; i++)
            b = b + a[i];    
    b = b/((T) n);
    return b;
}

template <typename T> void opuArraySetValue(T *y, T a, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

template <typename T> void opuArrayMultiplyScalar(T *y, T a, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = a*y[i];        
}

template <typename T> void opuArrayAddScalar(T *y, T a, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = y[i] + a;        
}

template <typename T> void opuArrayCopy(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = x[i];        
}

template <typename T> void opuArrayMinus(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = -x[i];        
}

template <typename T> void opuArrayAbs(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = fabs(x[i]);        
}

template <typename T> void opuArraySqrt(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = sqrt(x[i]);        
}

template <typename T> void opuArraySin(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = sin(x[i]);        
}

template <typename T> void opuArrayCos(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = cos(x[i]);        
}

template <typename T> void opuArrayTan(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = tan(x[i]);        
}

template <typename T> void opuArrayAsin(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = asin(x[i]);        
}

template <typename T> void opuArrayAcos(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = acos(x[i]);        
}

template <typename T> void opuArrayAtan(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = atan(x[i]);        
}

template <typename T> void opuArraySinh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = sinh(x[i]);        
}

template <typename T> void opuArrayCosh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = cosh(x[i]);        
}

template <typename T> void opuArrayTanh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = tanh(x[i]);        
}

template <typename T> void opuArrayAsinh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = asinh(x[i]);        
}

template <typename T> void opuArrayAcosh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = acosh(x[i]);        
}

template <typename T> void opuArrayAtanh(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = atanh(x[i]);        
}

template <typename T> void opuArrayExp(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = exp(x[i]);        
}

template <typename T> void opuArrayLog(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = log(x[i]);        
}

template <typename T> void opuArrayCeil(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = ceil(x[i]);        
}

template <typename T> void opuArrayFloor(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = floor(x[i]);        
}

template <typename T> void opuArrayErf(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = erf(x[i]);        
}

template <typename T> void opuArrayErfc(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = erfc(x[i]);        
}

template <typename T> void opuArraySquare(T *y, T *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = x[i]*x[i];        
}

template <typename T> void opuArrayPower(T *y, T *x, int p, int n)
{    
    for (int i=0; i<n; i++) {
        y[i] = x[i];    
        for (int j=1; j<p; j++)
            y[i] = y[i]*x[i];
    }
}

template <typename T> void opuArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    for (int i=0; i<n; i++)    
        C[i+n*i] = a*C[i+n*i];                
}

template <typename T> void opuArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{            
    for (int i=0; i<n; i++)    
        C[i+n*i] += a*x[i];                
}

template <typename T> void opuArrayRowAverage(T *y, T *x, int m, int n)
{    
    for (int j=0; j<n; j++)
    {        
        T avg = 0;
        int i;
        for (i=0; i<m; i++)
            avg = avg + x[i + m*j];
        avg = avg/((T) m);
        for (i=0; i<m; i++)
            y[i + m*j] = avg;         
    }        
}

template <typename T> void opuArrayAXPB(T *y, T *x, T a, T b, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = a*x[i]+b;        
}

template <typename T> void opuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{    
    for (int i=0; i<n; i++) 
        z[i] = a*x[i] + b*y[i];        
}

template <typename T> void opuArrayAXY(T *s, T *x, T *y, T a, int n)
{    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i];        
}

template <typename T> void opuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i]*z[i];        
}

template <typename T> void opuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i]*y[i] + b*z[i];        
}

template <typename T> void opuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{    
    for (int i=0; i<n; i++) 
        s[i] = a*x[i] + b*y[i] + c*z[i];        
}

template <typename T> void opuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+Q*k];
    }            
}

template <typename T> void opuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+Q*k] = un[idx];        
    }            
}

template <typename T> void opuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        C[i+I*j+I*J*s] = 0.0;
        for (int k=0; k<K; k++)
            C[i+I*j+I*J*s] += A[i+I*k+I*K*s]*B[k+K*j+K*J*s];
    }            
}

template <typename T> void opuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    for (int idx=0; idx<N; idx++)
    {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        for (int k=0; k<K; k++)
            C[idx] += A[i+I*k+Q*s]*B[k+K*j+P*s];
    }            
}

template <typename T> void opuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{    
    for (int i=0; i<nent; i++)         
    {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++)
            ucg[i] += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = ucg[i]/((T) nelem);
    }
}

template <typename T> void opuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{    
    for (int i=0; i<nent; i++)         
    {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            for (int j=0; j<npe; j++)
                ucg[i] += udg[j+npe*e]; 
        }
        ucg[i] = ucg[i]/((T) (nelem*npe));
    }
}

// template <typename T> void opuArrayCG2DG(T *udg, T *ucg, int *cgelcon, int npe, int ne)
// {        
//     for (int i=0; i<ne; i++)         
//         for (int k=0; k<npe; k++)
//             ucg[k+npe*i] = vcg[cgelcon[k+npe*i]];       
// }

template void opuGetArrayAtIndex(double*, double*, int*, int);
template void opuPutArrayAtIndex(double*, double*, int*, int);
template void opuArrayPlusXAtIndex(double*, double*, int*, int);
template void opuArrayMinusXAtIndex(double*, double*, int*, int);
template void opuArrayAXPYAtIndex(double*, double*, double, int*, int);
template void opuArrayAverageAtIndex(double*, double*, int*, int*, int);

template double opuArrayMin(double*, int);
template double opuArrayMax(double*, int);
template double opuArraySum(double*, int);
template double opuArraySquareSum(double*, int);
template double opuArrayMean(double*, int);
template void opuArraySetValue(double*, double, int);
template void opuArrayAddScalar(double*, double, int);
template void opuArrayMultiplyScalar(double*, double, int);

template void opuArrayCopy(double*, double*, int);
template void opuArrayMinus(double*, double*, int);
template void opuArrayAbs(double*, double*, int);
template void opuArraySqrt(double*, double*, int);
template void opuArraySin(double*, double*, int);
template void opuArrayCos(double*, double*, int);
template void opuArrayTan(double*, double*, int);
template void opuArrayAsin(double*, double*, int);
template void opuArrayAcos(double*, double*, int);
template void opuArrayAtan(double*, double*, int);
template void opuArraySinh(double*, double*, int);
template void opuArrayCosh(double*, double*, int);
template void opuArrayTanh(double*, double*, int);
template void opuArrayAsinh(double*, double*, int);
template void opuArrayAcosh(double*, double*, int);
template void opuArrayAtanh(double*, double*, int);
template void opuArrayExp(double*, double*, int);
template void opuArrayLog(double*, double*, int);
template void opuArrayCeil(double*, double*, int);
template void opuArrayFloor(double*, double*, int);
template void opuArrayErf(double*, double*, int);
template void opuArrayErfc(double*, double*, int);
template void opuArraySquare(double*, double*, int);
template void opuArrayPower(double*, double*, int, int);

template void opuArrayMultiplyScalarDiagonal(double*, double, int);
template void opuArrayAddVectorToDiagonal(double*, double*, double, int);
template void opuArrayRowAverage(double*, double*, int, int);
template void opuArrayAXPB(double*, double*, double, double, int);
template void opuArrayAXPBY(double*, double*, double*, double, double, int);
template void opuArrayAXY(double*, double*, double*, double, int);
template void opuArrayAXYZ(double*, double*, double*, double*, double, int);
template void opuArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void opuArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void opuArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void opuArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void opuArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void opuArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void opuArrayDG2CG(double*, double*, int*, int*, int);
template void opuArrayDG2CG2(double*, double*, int*, int*, int, int);

template void opuGetArrayAtIndex(float*, float*, int*, int);
template void opuPutArrayAtIndex(float*, float*, int*, int);
template void opuArrayPlusXAtIndex(float*, float*, int*, int);
template void opuArrayMinusXAtIndex(float*, float*, int*, int);
template void opuArrayAXPYAtIndex(float*, float*, float, int*, int);
template void opuArrayAverageAtIndex(float*, float*, int*, int*, int);

template float opuArrayMin(float*, int);
template float opuArrayMax(float*, int);
template float opuArraySum(float*, int);
template float opuArraySquareSum(float*, int);
template float opuArrayMean(float*, int);
template void opuArraySetValue(float*, float, int);
template void opuArrayAddScalar(float*, float, int);
template void opuArrayMultiplyScalar(float*, float, int);

template void opuArrayCopy(float*, float*, int);
template void opuArrayMinus(float*, float*, int);
template void opuArrayAbs(float*, float*, int);
template void opuArraySqrt(float*, float*, int);
template void opuArraySin(float*, float*, int);
template void opuArrayCos(float*, float*, int);
template void opuArrayTan(float*, float*, int);
template void opuArrayAsin(float*, float*, int);
template void opuArrayAcos(float*, float*, int);
template void opuArrayAtan(float*, float*, int);
template void opuArraySinh(float*, float*, int);
template void opuArrayCosh(float*, float*, int);
template void opuArrayTanh(float*, float*, int);
template void opuArrayAsinh(float*, float*, int);
template void opuArrayAcosh(float*, float*, int);
template void opuArrayAtanh(float*, float*, int);
template void opuArrayExp(float*, float*, int);
template void opuArrayLog(float*, float*, int);
template void opuArrayCeil(float*, float*, int);
template void opuArrayFloor(float*, float*, int);
template void opuArrayErf(float*, float*, int);
template void opuArrayErfc(float*, float*, int);
template void opuArraySquare(float*, float*, int);
template void opuArrayPower(float*, float*, int, int);

template void opuArrayMultiplyScalarDiagonal(float*, float, int);
template void opuArrayAddVectorToDiagonal(float*, float*, float, int);
template void opuArrayRowAverage(float*, float*, int, int);
template void opuArrayAXPB(float*, float*, float, float, int);
template void opuArrayAXPBY(float*, float*, float*, float, float, int);
template void opuArrayAXY(float*, float*, float*, float, int);
template void opuArrayAXYZ(float*, float*, float*, float*, float, int);
template void opuArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void opuArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void opuArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void opuArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void opuArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void opuArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void opuArrayDG2CG(float*, float*, int*, int*, int);
template void opuArrayDG2CG2(float*, float*, int*, int*, int, int);

template int opuArrayMin(int*, int);
template int opuArrayMax(int*, int);
template int opuArraySum(int*, int);
template int opuArraySquareSum(int*, int);
template void opuArrayCopy(int*, int*, int);
template void opuArraySetValue(int*, int, int);
template void opuArrayAddScalar(int*, int, int);
template void opuArrayMultiplyScalar(int*, int, int);
template void opuArrayAXPB(int*, int*, int, int, int);
template void opuArrayAXPBY(int*, int*, int*, int, int, int);
template void opuArrayExtract(int*, int*, int, int, int, int, int, int, int, int, int);
template void opuArrayInsert(int*, int*, int, int, int, int, int, int, int, int, int);

#endif


