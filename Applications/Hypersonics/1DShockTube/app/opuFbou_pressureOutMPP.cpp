template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T tau1 = tau[0];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T wdg6 = wdg[5*ng+i];
		T nlg1 = nlg[0*ng+i];
		T t2 = param1+param2+param3+param4+param5;
		T t3 = 1.0/t2;
		f[0*ng+i] = tau1*(param1-uhg1)+nlg1*param1*param6*t3;
		f[1*ng+i] = tau1*(param2-uhg2)+nlg1*param2*param6*t3;
		f[2*ng+i] = tau1*(param3-uhg3)+nlg1*param3*param6*t3;
		f[3*ng+i] = tau1*(param4-uhg4)+nlg1*param4*param6*t3;
		f[4*ng+i] = tau1*(param5-uhg5)+nlg1*param5*param6*t3;
		f[5*ng+i] = nlg1*(wdg6+(param6*param6)*t3)+tau1*(param6-uhg6);
		f[6*ng+i] = tau1*(param7-uhg7)+nlg1*param6*(param7*t3+t3*wdg6);
	}
}

template <typename T> void opuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T tau2 = tau[1];
		T tau3 = tau[2];
		T tau4 = tau[3];
		T tau5 = tau[4];
		T tau6 = tau[5];
		T tau7 = tau[6];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];

		int nspecies = 5;
		T rho_inf = 0.0085;
		T u_inf = 542.0;
		T rhoe_inf = rho_inf * u_inf * u_inf;
		double rhovec[5] = {udg1*rho_inf, udg2*rho_inf, udg3*rho_inf, udg4*rho_inf, udg5*rho_inf};

		// double rhovec[5] = {udg1, udg2, udg3, udg4, udg5};
		double rhoe = abs(udg7-(udg6*udg6)*(udg1/2.0+udg2/2.0+udg3/2.0+udg4/2.0+udg5/2.0)*1.0/pow(udg1+udg2+udg3+udg4+udg5,2.0));
		// std::cout << rhoe << std::endl;
		// if (rhoe < 0){
		// 	std::cout << xdg1 << std::endl;
		// }
		rhoe = rhoe * (rhoe_inf);
		mix->setState(rhovec, &rhoe, 0);

		T wdg6 = mix->P()/rhoe_inf;
		// T wdg6 = wdg[5*ng+i];
		T nlg1 = nlg[0*ng+i];
		T t2 = udg1+udg2+udg3+udg4+udg5;
		T t3 = 1.0/t2;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t3*udg1*udg6;
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t3*udg2*udg6;
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t3*udg3*udg6;
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t3*udg4*udg6;
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t3*udg5*udg6;
		f[5*ng+i] = nlg1*(wdg6+t3*(udg6*udg6))+tau6*(udg6-uhg6);
		f[6*ng+i] = tau7*(udg7-uhg7)+nlg1*udg6*(t3*udg7+t3*wdg6);
	}
}

template <typename T> void opuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
