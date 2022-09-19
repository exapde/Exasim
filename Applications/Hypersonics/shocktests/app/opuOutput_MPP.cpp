template <typename T> void opuOutput(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	// std::cout << "begin output" << std::endl;
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
		T udg8 = udg[j+npe*7+npe*nc*k];
		T udg9 = udg[j+npe*8+npe*nc*k];
		T udg10 = udg[j+npe*9+npe*nc*k];
		T udg11 = udg[j+npe*10+npe*nc*k];
		T udg12 = udg[j+npe*11+npe*nc*k];
		T udg13 = udg[j+npe*12+npe*nc*k];
		T udg14 = udg[j+npe*13+npe*nc*k];
		T odg1 = odg[j+npe*0+npe*nco*k];
		T odg2 = odg[j+npe*1+npe*nco*k];

		T t1pi = 1.0/3.141592653589793;
		udg1 = udg1*(t1pi*atan(udg1*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg2 = udg2*(t1pi*atan(udg2*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg3 = udg3*(t1pi*atan(udg3*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg4 = udg4*(t1pi*atan(udg4*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg5 = udg5*(t1pi*atan(udg5*1.0E+3)+1.0/2.0)+3.183097800805168E-4;


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
		T uinf1 = mix->P();

		f[j+npe*0+npe*nce*k] = uinf1;
		f[j+npe*1+npe*nce*k] = odg1;
	}
	// std::cout << "end output" << std::endl;
}

template void opuOutput(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuOutput(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
