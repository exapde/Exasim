template <typename T> void opuOutput(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T udg1 = udg[j+npe*0+npe*nc*k];
		T udg2 = udg[j+npe*1+npe*nc*k];
		T udg3 = udg[j+npe*2+npe*nc*k];
		T udg4 = udg[j+npe*3+npe*nc*k];
		T udg5 = udg[j+npe*4+npe*nc*k];
		T udg6 = udg[j+npe*5+npe*nc*k];
		T udg7 = udg[j+npe*6+npe*nc*k];

        int nspecies = 5;
		double rho_scale = uinf[0];
		double u_scale = uinf[1];
		double rhoe_scale = uinf[2];

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);

		f[j+npe*0+npe*nce*k] = mix->P();
		f[j+npe*1+npe*nce*k] = mix->T();
	}
}

template void opuOutput(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuOutput(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
