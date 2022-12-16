template <typename T> void opuOutput(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
{
    double Uwork[5];
    double dTdU[6];
    double D_i[5];
    double h_i[5];
    int nspecies = 5;
    int ndim = 1;
    double rho_scale = uinf[0];
    double u_scale = uinf[1];
    double rhoe_scale = uinf[2];
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param15 = param[14];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf3 = uinf[2];
		T uinf4 = uinf[3];
		T uinf5 = uinf[4];
		T uinf6 = uinf[5];
		T uinf7 = uinf[6];
		T uinf8 = uinf[7];
		T uinf9 = uinf[8];
		T uinf10 = uinf[9];
		T uinf11 = uinf[10];
		T uinf12 = uinf[11];
		T uinf13 = uinf[12];
		T uinf14 = uinf[13];
		T uinf15 = uinf[14];
		T uinf16 = uinf[15];
		T uinf17 = uinf[16];
		T uinf18 = uinf[17];
		T uinf19 = uinf[18];
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
        udg1 = udg1*(t1pi*atan(udg1*1.0E+12)+1.0/2.0)+3.182454300088011e-13;
		udg2 = udg2*(t1pi*atan(udg2*1.0E+12)+1.0/2.0)+3.182454300088011e-13;
		udg3 = udg3*(t1pi*atan(udg3*1.0E+12)+1.0/2.0)+3.182454300088011e-13;
		udg4 = udg4*(t1pi*atan(udg4*1.0E+12)+1.0/2.0)+3.182454300088011e-13;
		udg5 = udg5*(t1pi*atan(udg5*1.0E+12)+1.0/2.0)+3.182454300088011e-13;

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies, 1);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);
		uinf1 = mix->P()/rhoe_scale;
        uinf2 = mix->mixtureFrozenGamma();
        dT_dUstate(dTdU, Ustate, Uwork, nspecies, ndim, mix);
        nondimensionalize_dT_dUstate(dTdU, (double*)uinf, nspecies, ndim);
        uinf13 = dTdU[0]; // dT_drho_1
		uinf14 = dTdU[1]; // dT_drho_2
		uinf15 = dTdU[2]; // dT_drho_3
		uinf16 = dTdU[3]; // dT_drho_4
		uinf17 = dTdU[4]; // dT_drho_5
		uinf18 = dTdU[5]; // dT_drhoe

        mix->averageDiffusionCoeffs(D_i);
        nondimensionalize_diffusionCoeffs(D_i, (double*)uinf, nspecies, ndim);
        mix->speciesHOverRT(h_i);
        nondimensionalize_enthalpies(h_i, (double*)uinf, nspecies, ndim);

		uinf2  = D_i[0]; //  D_1
		uinf3  = D_i[1]; //  D_2
		uinf4  = D_i[2]; //  D_3 
		uinf5  = D_i[3]; //  D_4
		uinf6  = D_i[4]; //  D_5
		uinf7  = h_i[0]; //  h_1
		uinf8  = h_i[1]; //  h_2
		uinf9  = h_i[2]; //  h_3
		uinf10 = h_i[3]; //  h_4
		uinf11 = h_i[4]; //  h_5 
		uinf12 = mix->viscosity() / uinf[4]; //  mu
		uinf13 = mix->frozenThermalConductivity() / uinf[5]; //  kappa

        
		f[j+npe*0+npe*nce*k] = uinf2;
		f[j+npe*1+npe*nce*k] = uinf3;
		f[j+npe*2+npe*nce*k] = uinf4;
		f[j+npe*3+npe*nce*k] = uinf5;
		f[j+npe*4+npe*nce*k] = uinf6;
		f[j+npe*5+npe*nce*k] = uinf7;
		f[j+npe*6+npe*nce*k] = uinf8;
		f[j+npe*7+npe*nce*k] = uinf9;
		f[j+npe*8+npe*nce*k] = uinf10;
		f[j+npe*9+npe*nce*k] = uinf11;
        f[j+npe*10+npe*nce*k] = uinf12;
        f[j+npe*11+npe*nce*k] = uinf13;
	}
}

template void opuOutput(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuOutput(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
