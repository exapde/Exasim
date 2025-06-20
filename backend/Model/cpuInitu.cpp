void cpuInitu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)
{

  for (int i = 0; i < N; ++i) {
    int j = i%npe; 
    int k = i/npe; 
    dstype mu4 = mu[4];
    dstype mu5 = mu[5];
    dstype mu6 = mu[6];
    dstype mu7 = mu[7];


    f[j+npe*0 +npe*nce*k] = mu4;
    f[j+npe*1 +npe*nce*k] = mu5;
    f[j+npe*2 +npe*nce*k] = mu6;
    f[j+npe*3 +npe*nce*k] = mu7;
  }
}

