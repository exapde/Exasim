#ifndef __CPUCORE_H__
#define __CPUCORE_H__

template <typename T> void cpuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2);
void cpuIndexPermute12(int *index, int I1, int I2, int I3);
void cpuIndexPermute13(int *index, int I1, int I2, int I3, int I4);
void cpuIndexPermute23(int *index, int I1, int I2, int I3, int I4);
template <typename T> void cpuPermute(T *B, T *A, int *index, int N);
template <typename T> void cpuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM);
template <typename T> void cpuPermute12(T *B, T *A, int I1, int I2, int I3);
template <typename T> void cpuPermute13(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void cpuPermute23(T *B, T *A, int I1, int I2, int I3, int I4);

template <typename T> void cpuGetArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuPutArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n);

template <typename T> T cpuArrayMin(T *a, int n);
template <typename T> T cpuArrayMax(T *a, int n);
template <typename T> T cpuArraySum(T *a, int n);
template <typename T> T cpuArraySquareSum(T *a, int n);
template <typename T> T cpuArrayMean(T *a, int n);
template <typename T> void cpuArraySetValue(T *y, T a, int n);
template <typename T> void cpuArrayAddScalar(T *y, T a, int n);
template <typename T> void cpuArrayMultiplyScalar(T *y, T a, int n);

template <typename T> void cpuArrayCopy(T *y, T *x, int n);
template <typename T> void cpuArrayMinus(T *y, T *x, int n);
template <typename T> void cpuArrayAbs(T *y, T *x, int n);
template <typename T> void cpuArraySqrt(T *y, T *x, int n);
template <typename T> void cpuArraySin(T *y, T *x, int n);
template <typename T> void cpuArrayCos(T *y, T *x, int n);
template <typename T> void cpuArrayTan(T *y, T *x, int n);
template <typename T> void cpuArrayAsin(T *y, T *x, int n);
template <typename T> void cpuArrayAcos(T *y, T *x, int n);
template <typename T> void cpuArrayAtan(T *y, T *x, int n);
template <typename T> void cpuArraySinh(T *y, T *x, int n);
template <typename T> void cpuArrayCosh(T *y, T *x, int n);
template <typename T> void cpuArrayTanh(T *y, T *x, int n);
template <typename T> void cpuArrayAsinh(T *y, T *x, int n);
template <typename T> void cpuArrayAcosh(T *y, T *x, int n);
template <typename T> void cpuArrayAtanh(T *y, T *x, int n);
template <typename T> void cpuArrayExp(T *y, T *x, int n);
template <typename T> void cpuArrayLog(T *y, T *x, int n);
template <typename T> void cpuArrayCeil(T *y, T *x, int n);
template <typename T> void cpuArrayFloor(T *y, T *x, int n);
template <typename T> void cpuArrayErf(T *y, T *x, int n);
template <typename T> void cpuArrayErfc(T *y, T *x, int n);
template <typename T> void cpuArraySquare(T *y, T *x, int n);
template <typename T> void cpuArrayPower(T *y, T *x, int p, int n);

template <typename T> void cpuArrayMultiplyScalarDiagonal(T *C, T a, int n);
template <typename T> void cpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n);
template <typename T> void cpuArrayRowAverage(T *a, T *b, int m, int n);
template <typename T> void cpuArrayAXPB(T *y, T *x, T a, T b, int n);
template <typename T> void cpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n);
template <typename T> void cpuArrayAXY(T *s, T *x, T *y, T a, int n);
template <typename T> void cpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n);
template <typename T> void cpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n);
template <typename T> void cpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n);
template <typename T> void cpuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void cpuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void cpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void cpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void cpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent);
template <typename T> void cpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe);

template <typename T> void cpuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts);

template <typename T> void cpuElemGeom(T *Xx, T *jac, T *Jg, int ne, int ng, int nd);
template <typename T> void cpuFaceGeom(T *nlg, T *jacg, T *Jg, int nf, int ng, int nd);
template <typename T> void cpuElemGeom1D(T *jac, T *Xx, T *Jg, int nga);
template <typename T> void cpuElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22, T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga);
template <typename T> void cpuElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33, T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga);
template <typename T> void cpuFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga);
template <typename T> void cpuFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga);
template <typename T> void cpuFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga);

template <typename T> void cpuApplyJac(T *sg, T *fhg, T *jac, int nga, int ncu, int ngf);
template <typename T> void cpuApplyJac1(T *sg, T *jac, int nga, int ncu);
template <typename T> void cpuApplyJac2(T *sg, T *jac, T *ug, T *su, T *fc_u, int nga, int ncu);
template <typename T> void cpuApplyXx1(T *sg, T *ug, T *Xx, int nga, int nd, int ncu);
template <typename T> void cpuApplyXx2(T *sg, T *fg, T *Xx, int nga, int nd, int ncu);
template <typename T> void cpuApplyXx3(T *sg, T *ug, T *Xx, int nge, int nd, int ncu, int ne);
template <typename T> void cpuApplyXx4(T *rg, T *sg, T *fg, T *Xx, T *jac, int nge, int nd, int ncu, int ne);
template <typename T> void cpuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd);
template <typename T> void cpuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd, int ngf);
template <typename T> void cpuApplyFactor(T *Rfac, T *R, T *fac, int npe, int M, int N);
template <typename T> void cpuApplyFactorJac(T *Rfac, T *R, T *fac, T *jac, int npe, int M, int N);
template <typename T> void cpuApplyJacInv(T *Rfac, T *R, T *jac, int M, int N);
template <typename T> void cpuApplyJac(T *Rfac, T *R, T *jac, int M, int N);
template <typename T> void cpuShapJac(T *shapjac, T *shapegt, T *jac, int nge, int M, int N);

template <typename T> void cpuAverageFlux(T *fg, int N);
template <typename T> void cpuAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd);
template <typename T> void cpuAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M);
template <typename T> void cpuAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints);
template <typename T> void cpuAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu);

#endif  

