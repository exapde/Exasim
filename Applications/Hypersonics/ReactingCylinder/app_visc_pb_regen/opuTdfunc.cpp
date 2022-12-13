template <typename T> void opuTdfunc(T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		f[0*ng+i] = 1.0;
		f[1*ng+i] = 1.0;
		f[2*ng+i] = 1.0;
		f[3*ng+i] = 1.0;
		f[4*ng+i] = 1.0;
		f[5*ng+i] = 1.0;
		f[6*ng+i] = 1.0;
		f[7*ng+i] = 1.0;
	}
}

template void opuTdfunc(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuTdfunc(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
