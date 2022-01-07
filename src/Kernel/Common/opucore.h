#ifndef __OPUCORE_H__
#define __OPUCORE_H__

template <typename T> void opuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2);
void opuIndexPermute12(int *index, int I1, int I2, int I3);
void opuIndexPermute13(int *index, int I1, int I2, int I3, int I4);
void opuIndexPermute23(int *index, int I1, int I2, int I3, int I4);
template <typename T> void opuPermute(T *B, T *A, int *index, int N);
template <typename T> void opuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM);
template <typename T> void opuPermute12(T *B, T *A, int I1, int I2, int I3);
template <typename T> void opuPermute13(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void opuPermute23(T *B, T *A, int I1, int I2, int I3, int I4);

template <typename T> void opuGetArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuPutArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayPlusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayMinusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n);
template <typename T> void opuArrayAverageAtIndex(T *y, T *x, int *ind1, int *ind2, int n);

template <typename T> T opuArrayMin(T *a, int n);
template <typename T> T opuArrayMax(T *a, int n);
template <typename T> T opuArraySum(T *a, int n);
template <typename T> T opuArraySquareSum(T *a, int n);
template <typename T> T opuArrayMean(T *a, int n);
template <typename T> void opuArraySetValue(T *y, T a, int n);
template <typename T> void opuArrayAddScalar(T *y, T a, int n);
template <typename T> void opuArrayMultiplyScalar(T *y, T a, int n);

template <typename T> void opuArrayCopy(T *y, T *x, int n);
template <typename T> void opuArrayMinus(T *y, T *x, int n);
template <typename T> void opuArrayAbs(T *y, T *x, int n);
template <typename T> void opuArraySqrt(T *y, T *x, int n);
template <typename T> void opuArraySin(T *y, T *x, int n);
template <typename T> void opuArrayCos(T *y, T *x, int n);
template <typename T> void opuArrayTan(T *y, T *x, int n);
template <typename T> void opuArrayAsin(T *y, T *x, int n);
template <typename T> void opuArrayAcos(T *y, T *x, int n);
template <typename T> void opuArrayAtan(T *y, T *x, int n);
template <typename T> void opuArraySinh(T *y, T *x, int n);
template <typename T> void opuArrayCosh(T *y, T *x, int n);
template <typename T> void opuArrayTanh(T *y, T *x, int n);
template <typename T> void opuArrayAsinh(T *y, T *x, int n);
template <typename T> void opuArrayAcosh(T *y, T *x, int n);
template <typename T> void opuArrayAtanh(T *y, T *x, int n);
template <typename T> void opuArrayExp(T *y, T *x, int n);
template <typename T> void opuArrayLog(T *y, T *x, int n);
template <typename T> void opuArrayCeil(T *y, T *x, int n);
template <typename T> void opuArrayFloor(T *y, T *x, int n);
template <typename T> void opuArrayErf(T *y, T *x, int n);
template <typename T> void opuArrayErfc(T *y, T *x, int n);
template <typename T> void opuArraySquare(T *y, T *x, int n);
template <typename T> void opuArrayPower(T *y, T *x, int p, int n);

template <typename T> void opuArrayMultiplyScalarDiagonal(T *C, T a, int n);
template <typename T> void opuArrayAddVectorToDiagonal(T *C, T *x, T a, int n);
template <typename T> void opuArrayRowAverage(T *a, T *b, int M, int N);
template <typename T> void opuArrayAXPB(T *y, T *x, T a, T b, int n);
template <typename T> void opuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n);
template <typename T> void opuArrayAXY(T *s, T *x, T *y, T a, int n);
template <typename T> void opuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n);
template <typename T> void opuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n);
template <typename T> void opuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n);
template <typename T> void opuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void opuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void opuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void opuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void opuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent);
template <typename T> void opuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe);

template <typename T> void opuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts);

template <typename T> void opuElemGeom(T *Xx, T *jac, T *Jg, int ne, int ng, int nd);
template <typename T> void opuFaceGeom(T *nlg, T *jacg, T *Jg, int nf, int ng, int nd);
template <typename T> void opuElemGeom1D(T *jac, T *Xx, T *Jg, int nga);
template <typename T> void opuElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22, T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga);
template <typename T> void opuElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33, T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga);
template <typename T> void opuFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga);
template <typename T> void opuFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga);
template <typename T> void opuFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga);
template <typename T> void opuFixNormal1D(T *nlg, int *facecon, int na);

template <typename T> void opuApplyJac(T *sg, T *fhg, T *jac, int nga, int ncu, int ngf);
template <typename T> void opuApplyJac1(T *sg, T *jac, int nga, int ncu);
template <typename T> void opuApplyJac2(T *sg, T *jac, T *ug, T *su, T *fc_u, int nga, int ncu);
template <typename T> void opuApplyXx1(T *sg, T *ug, T *Xx, int nga, int nd, int ncu);
template <typename T> void opuApplyXx2(T *sg, T *fg, T *Xx, int nga, int nd, int ncu);
template <typename T> void opuApplyXx3(T *sg, T *ug, T *Xx, int nge, int nd, int ncu, int ne);
template <typename T> void opuApplyXx4(T *rg, T *sg, T *fg, T *Xx, T *jac, int nge, int nd, int ncu, int ne);
template <typename T> void opuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd);
template <typename T> void opuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd, int ngf);
template <typename T> void opuApplyFactor(T *Rfac, T *R, T *fac, int npe, int M, int N);
template <typename T> void opuApplyFactorJac(T *Rfac, T *R, T *fac, T *jac, int npe, int M, int N);
template <typename T> void opuApplyJacInv(T *Rfac, T *R, T *jac, int M, int N);
template <typename T> void opuApplyJac(T *Rfac, T *R, T *jac, int M, int N);
template <typename T> void opuShapJac(T *shapjac, T *shapegt, T *jac, int nge, int M, int N);

template <typename T> void opuAverageFlux(T *fg, int N);
template <typename T> void opuAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd);
template <typename T> void opuAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M);
template <typename T> void opuAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints);
template <typename T> void opuAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu);

template <typename T> void opuStgHomoTurb(T *up, T *xdg, T *stgdata, T *uc, T t, int M, int N, int nd);
template <typename T> void opuStgHomoTurb(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N, int nd);
template <typename T> void opuStgHomoTurb2(T *up, T *xdg, T *stgdata, T *ui, T *uc, T *c, T *param, T t, int M, int N, int nd);

#endif  

