template <typename T> void opuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		// T xdg1 = xdg[j+npe*0+npe*ncx*k];
		// T udg1, udg2, udg3, udg4, udg5, udg6, udg7;
		// T param1 = param[0];
		// T param2 = param[1];
		// T param3 = param[2];
		// T param4 = param[3];
		// T param5 = param[4];
		// T param6 = param[5];
		// T param7 = param[6];
		// T param8 = param[7];
		// T param9 = param[8];
		// T param10 = param[9];
		// T param11 = param[10];
		// T param12 = param[11];
		// T param13 = param[12];
		// T param14 = param[13];
		// // if (xdg1 < 0)
        // // {
		///////// Inflow version
		// udg1 = param1;
		// udg2 = param2;
		// udg3 = param3;
		// udg4 = param4;
		// udg5 = param5;
		// udg6 = param6;
		// udg7 = param7;
        // }
		///////// Conditional for full domain
        // else
        // {
        //     udg1 = param8;
        //     udg2 = param9;
        //     udg3 = param10;
        //     udg4 = param11;
        //     udg5 = param12;
        //     udg6 = param13;
		// 	udg7 = param14;
        // }

		///////// Smooth out initial
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T t2 = param1/2.0;
		T t3 = param2/2.0;
		T t4 = param3/2.0;
		T t5 = param4/2.0;
		T t6 = param5/2.0;
		T t7 = param6/2.0;
		T t8 = param7/2.0;
		T t9 = param8/2.0;
		T t10 = param9/2.0;
		T t11 = param10/2.0;
		T t12 = param11/2.0;
		T t13 = param12/2.0;
		T t14 = param13/2.0;
		T t15 = param14/2.0;
		T t16 = xdg1*1.0E+3;
		T t17 = tanh(t16);
		T udg1 = t2+t9+t17*(t2-t9);
		T udg2 = t3+t10+t17*(t3-t10);
		T udg3 = t4+t11+t17*(t4-t11);
		T udg4 = t5+t12+t17*(t5-t12);
		T udg5 = t6+t13+t17*(t6-t13);
		T udg6 = t7+t14+t17*(t7-t14);
		T udg7 = t8+t15+t17*(t8-t15);
		int nspecies = 5;
		T rho_inf = 0.0085;
		T u_inf = 542.0;
		T rhoe_inf = rho_inf * u_inf * u_inf;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};

		// double rhovec[5] = {udg1, udg2, udg3, udg4, udg5};
		double wdot[5];
		double rhoe = udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0);
		// std::cout << rhoe << std::endl;
		rhoe = rhoe * rhoe_inf;
		mix->setState(rhovec, &rhoe, 0);
		mix->netProductionRates(wdot);
		for (int ispecies = 0; ispecies < nspecies; ispecies++ )
		{
			f[j+npe*ispecies+npe*nce*k] = wdot[ispecies];
		}
		f[j+npe*nspecies+npe*nce*k] = mix->P();

	}
}

template void opuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int, Mutation::Mixture *);
template void opuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int, Mutation::Mixture *);
