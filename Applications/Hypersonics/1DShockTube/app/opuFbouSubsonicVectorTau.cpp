template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		// Try longer time
		// Try fixing inlet pressure
		// Using exact solution 
		// Average inlet flux
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
		T nlg1 = nlg[0*ng+i];

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
        T uinf4 = mix->P()/rhoe_scale;
		// T uinf4 = 16894.0/rhoe_scale;

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
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t29*udg6*(t23+3.183097800805168E-4);
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t29*udg6*(t24+3.183097800805168E-4);
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t29*udg6*(t25+3.183097800805168E-4);
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t29*udg6*(t26+3.183097800805168E-4);
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t29*udg6*(t27+3.183097800805168E-4);
		f[5*ng+i] = nlg1*(uinf4+t29*(udg6*udg6))+tau6*(udg6-uhg6);
		f[6*ng+i] = tau7*(udg7-uhg7)+nlg1*udg6*(t29*udg7+t29*uinf4);
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
		T nlg1 = nlg[0*ng+i];

        int nspecies = 5;
		double rho_scale = uinf[0];
		double u_scale = uinf[1];
		double rhoe_scale = uinf[2];

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
        // double Ucons[7] = {uhg1, uhg2, uhg3, uhg4, uhg5, uhg6, uhg7};
		double Ustate[6];
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);
        // T uinf4 = param[3*nspecies+6]/rhoe_scale;
        T uinf4 = mix->P()/rhoe_scale;

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

		T fnudg1 = nlg1*t29*udg6*(t23+3.183097800805168E-4);
		T fnudg2 = nlg1*t29*udg6*(t24+3.183097800805168E-4);
		T fnudg3 = nlg1*t29*udg6*(t25+3.183097800805168E-4);
		T fnudg4 = nlg1*t29*udg6*(t26+3.183097800805168E-4);
		T fnudg5 = nlg1*t29*udg6*(t27+3.183097800805168E-4);
		T fnudg6 = nlg1*(uinf4+t29*(udg6*udg6));
		T fnudg7 = nlg1*udg6*(t29*udg7+t29*uinf4);

		T t3h = uhg1*1.0E+3;
		T t4h = uhg2*1.0E+3;
		T t5h = uhg3*1.0E+3;
		T t6h = uhg4*1.0E+3;
		T t7h = uhg5*1.0E+3;
		T t8h = atan(t3h);
		T t9h = atan(t4h);
		T t10h = atan(t5h);
		T t11h = atan(t6h);
		T t12h = atan(t7h);
		T t13h = t2*t8h;
		T t14h = t2*t9h;
		T t15h = t2*t10h;
		T t16h = t2*t11h;
		T t17h = t2*t12h;
		T t18h = t13h+1.0/2.0;
		T t19h = t14h+1.0/2.0;
		T t20h = t15h+1.0/2.0;
		T t21h = t16h+1.0/2.0;
		T t22h = t17h+1.0/2.0;
		T t23h = t18h*uhg1;
		T t24h = t19h*uhg2;
		T t25h = t20h*uhg3;
		T t26h = t21h*uhg4;
		T t27h = t22h*uhg5;
		T t28h = t23h+t24h+t25h+t26h+t27h+1.591548900402584E-3;
		T t29h = 1.0/t28h;
		
		// T fnug1 = nlg1*t3*udg1*udg6;
		// T fnudg2  = nlg1*t3*udg2*udg6;
		// T fnudg3 = nlg1*t3*udg3*udg6;
		// T fnudg4 = nlg1*t3*udg4*udg6;
		// T fnudg5 = nlg1*t3*udg5*udg6;
		// T fnudg6 = nlg1*(wdg6+t3*(udg6*udg6));
		// T fnudg7 = nlg1*udg6*(t3*udg7+t3*wdg6);
		uinf4 = param[3*nspecies+6]/rhoe_scale;
		T fnuhg1 = nlg1*t29h*uhg6*(t23h+3.183097800805168E-4);
		T fnuhg2 = nlg1*t29h*uhg6*(t24h+3.183097800805168E-4);
		T fnuhg3 = nlg1*t29h*uhg6*(t25h+3.183097800805168E-4);
		T fnuhg4 = nlg1*t29h*uhg6*(t26h+3.183097800805168E-4);
		T fnuhg5 = nlg1*t29h*uhg6*(t27h+3.183097800805168E-4);
		T fnuhg6 = nlg1*(uinf4+t29h*(uhg6*uhg6));
		T fnuhg7 = nlg1*uhg6*(t29h*uhg7+t29h*uinf4);

		// f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t29*udg6*(t23+3.183097800805168E-4);
		// f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t29*udg6*(t24+3.183097800805168E-4);
		// f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t29*udg6*(t25+3.183097800805168E-4);
		// f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t29*udg6*(t26+3.183097800805168E-4);
		// f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t29*udg6*(t27+3.183097800805168E-4);
		// f[5*ng+i] = nlg1*(uinf4+t29*(udg6*udg6))+tau6*(udg6-uhg6);
		// f[6*ng+i] = tau7*(udg7-uhg7)+nlg1*udg6*(t29*udg7+t29*uinf4);

		f[0*ng+i] = tau1*(udg1-uhg1)+ 0.5 * (fnudg1 + fnuhg1);
		f[1*ng+i] = tau2*(udg2-uhg2)+0.5 * (fnudg2 + fnuhg2);
		f[2*ng+i] = tau3*(udg3-uhg3)+0.5 * (fnudg3 + fnuhg3);
		f[3*ng+i] = tau4*(udg4-uhg4)+0.5 * (fnudg4 + fnuhg4);
		f[4*ng+i] = tau5*(udg5-uhg5)+0.5 * (fnudg5 + fnuhg5);
		f[5*ng+i] = 0.5 * (fnudg6 + fnuhg6) +tau6*(udg6-uhg6);
		f[6*ng+i] = tau7*(udg7-uhg7)+0.5 * (fnudg7 + fnuhg7);
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
