void cpuInitu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)
{

  for (int i = 0; i < N; ++i) {
    int j = i%npe; 
    int k = i/npe; 


    f[j+npe*0 +npe*nce*k] = 0.0;
  }
}

