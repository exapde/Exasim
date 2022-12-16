template <typename T> void opuSourcew(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
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
		T wdg1 = wdg[j+npe*0+npe*ncw*k];
		T wdg2 = wdg[j+npe*1+npe*ncw*k];
		T wdg3 = wdg[j+npe*2+npe*ncw*k];
		T wdg4 = wdg[j+npe*3+npe*ncw*k];
		T wdg5 = wdg[j+npe*4+npe*ncw*k];
		T wdg6 = wdg[j+npe*5+npe*ncw*k];

		int nspecies = 5;
		T rho_inf = 0.0085;
		T u_inf = 542.0;
		T rhoe_inf = rho_inf * u_inf * u_inf;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};

		// double rhovec[5] = {udg1, udg2, udg3, udg4, udg5};
		double wdot[5];
		double rhoe = abs(udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0));
		// std::cout << rhoe << std::endl;
		// if (rhoe < 0){
			// std::cout << xdg1 << std::endl;
		// }
		rhoe = rhoe * rhoe_inf;
		rhoe = fmax(rhoe,200.0);
		// std::cout << rhoe << std::endl;
		mix->setState(rhovec, &rhoe, 0);
		mix->netProductionRates(wdot);
		for (int ispecies = 0; ispecies < nspecies; ispecies++ )
		{
			f[j+npe*ispecies+npe*nce*k] = wdot[ispecies]/(rho_inf*u_inf);
		}
		f[j+npe*nspecies+npe*nce*k] = mix->P()/rhoe_inf;
	}
}

template void opuSourcew(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuSourcew(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
