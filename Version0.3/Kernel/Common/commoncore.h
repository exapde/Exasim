#ifndef __COMMONCORE_H__
#define __COMMONCORE_H__

void faceperm(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf);
void faceperm2(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf);
void faceperm1(int *ind1, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf);
void elemperm(int *ind, int *indpts, int *eblks, int npe, int nc, int ncu, int nbe);

void faceperm(unsigned int *ind1, unsigned int *ind2, unsigned int *indpts, unsigned int *facecon, unsigned int *fblks, unsigned int npf, unsigned int ncu, unsigned int npe, unsigned int nc, unsigned int nbf);
void faceperm1(unsigned int *ind1, unsigned int *indpts, unsigned int *facecon, unsigned int *fblks, unsigned int npf, unsigned int ncu, unsigned int npe, unsigned int nc, unsigned int nbf);
void elemperm(unsigned int *ind, unsigned int *indpts, unsigned int *eblks, unsigned int npe, unsigned int nc, unsigned int ncu, unsigned int nbe);

template <typename T> T cpuArrayGetElementAtIndex(T *y, int n);
template <typename T> void cpuArraySetValueAtIndex(T *y, T a, int n);
template <typename T> void cpuApplyGivensRotation(T *H, T *s, T *cs, T *sn,  int i);
template <typename T> void cpuBackSolve(T *y, T *H, T *s, int i, int n);

#endif




