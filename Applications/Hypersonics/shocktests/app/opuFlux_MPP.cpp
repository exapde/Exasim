template <typename T> void opuFlux(T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
		// std::cout << "begin flux" << std::endl;

	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T udg13 = udg[12*ng+i];
		T udg14 = udg[13*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
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
		T uinf1 = mix->P()/rhoe_scale;

		T t2 = 1.0/3.141592653589793;
		T t3 = udg1*1.0E+3;
		T t4 = udg2*1.0E+3;
		T t5 = udg3*1.0E+3;
		T t6 = udg4*1.0E+3;
		T t7 = udg5*1.0E+3;
		T t8 = atan(t3);
		T t9 = atan(t4);
		T t10 = atan(t5);
		T t11 = atan(t6);
		T t12 = atan(t7);
		T t13 = t2*t8;
		T t14 = t2*t9;
		T t15 = t2*t10;
		T t16 = t2*t11;
		T t17 = t2*t12;
		T t18 = t13+1.0/2.0;
		T t19 = t14+1.0/2.0;
		T t20 = t15+1.0/2.0;
		T t21 = t16+1.0/2.0;
		T t22 = t17+1.0/2.0;
		T t23 = t18*udg1;
		T t24 = t19*udg2;
		T t25 = t20*udg3;
		T t26 = t21*udg4;
		T t27 = t22*udg5;
		T t28 = t23+t24+t25+t26+t27+1.591548900402584E-3;
		T t29 = 1.0/t28;
		f[0*ng+i] = odg1*udg8*(t18+(t2*t3)/((udg1*udg1)*1.0E+6+1.0))+t29*udg6*(t23+3.183097800805168E-4);
		f[1*ng+i] = odg1*udg9*(t19+(t2*t4)/((udg2*udg2)*1.0E+6+1.0))+t29*udg6*(t24+3.183097800805168E-4);
		f[2*ng+i] = odg1*udg10*(t20+(t2*t5)/((udg3*udg3)*1.0E+6+1.0))+t29*udg6*(t25+3.183097800805168E-4);
		f[3*ng+i] = odg1*udg11*(t21+(t2*t6)/((udg4*udg4)*1.0E+6+1.0))+t29*udg6*(t26+3.183097800805168E-4);
		f[4*ng+i] = odg1*udg12*(t22+(t2*t7)/((udg5*udg5)*1.0E+6+1.0))+t29*udg6*(t27+3.183097800805168E-4);
		f[5*ng+i] = uinf1+odg1*udg13+t29*(udg6*udg6);
		f[6*ng+i] = odg1*udg14+udg6*(t29*udg7+t29*uinf1);
	}
		// std::cout << "end flux" << std::endl;

}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
