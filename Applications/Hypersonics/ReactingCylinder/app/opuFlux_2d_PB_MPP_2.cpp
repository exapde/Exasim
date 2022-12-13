template <typename T> void opuFlux(T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
    double Ustate[6];
	int nspecies = 5;
	int ndim = 2;
	double rho_scale = uinf[0];
	double u_scale = uinf[1];
	double rhoe_scale = uinf[2];
	double omega_scale = rho_scale*u_scale;
    double mu_scale = uinf[4];
    double kappa_scale = uinf[5];
    double Uwork[5];
    double dTdU[6];
	double D_i[5];
    double h_i[5];
    
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param19 = param[18];
		T param20 = param[19];
		T param21 = param[20];
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
		T uinf20 = uinf[19];
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
		T udg15 = udg[14*ng+i];
		T udg16 = udg[15*ng+i];
		T udg17 = udg[16*ng+i];
		T udg18 = udg[17*ng+i];
		T udg19 = udg[18*ng+i];
		T udg20 = udg[19*ng+i];
		T udg21 = udg[20*ng+i];
		T udg22 = udg[21*ng+i];
		T udg23 = udg[22*ng+i];
		T udg24 = udg[23*ng+i];
		T odg1 = odg[0*ng+i];

		        T t1pi = 1.0/3.141592653589793;
        
        udg1 = udg1*(t1pi*atan(udg1*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg2 = udg2*(t1pi*atan(udg2*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg3 = udg3*(t1pi*atan(udg3*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg4 = udg4*(t1pi*atan(udg4*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg5 = udg5*(t1pi*atan(udg5*1.0E+12)+1.0/2.0)+3.182454300088011E-13;

		double Ucons[8] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7, udg8};
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, ndim);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies, ndim);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

        mix->setState(rhovec, &rhoe, 0);

		// mix->averageDiffusionCoeffs(D_i);
        // nondimensionalize_diffusionCoeffs(D_i, (double*)uinf, nspecies, ndim);
        mix->speciesHOverRT(h_i);
        nondimensionalize_enthalpies(h_i, (double*)uinf, nspecies, ndim);

        uinf1 = mix->P() / rhoe_scale;

        //D_i
        uinf2  = 0.0; uinf3  = 0.0; uinf4  = 0.0; uinf5  = 0.0; uinf6  = 0.0;
        // uinf2  = D_i[0]; uinf3  = D_i[1]; uinf4  = D_i[2]; uinf5  = D_i[3]; uinf6  = D_i[4]; 

        //h_i
        uinf7  = 0.0; uinf8  = 0.0; uinf9  = 0.0; uinf10 = 0.0; uinf11 = 0.0; 
		// uinf7  = h_i[0]; uinf8  = h_i[1]; uinf9  = h_i[2]; uinf10 = h_i[3]; uinf11 = h_i[4]; //  h_5 

		uinf20 = 1.0/0.9;
        uinf12 = 0.0; uinf13 = 0.0; //  kappa
		// uinf12 = mix->viscosity() / mu_scale; uinf13 = mix->frozenThermalConductivity() / kappa_scale; //  kappa

        dT_dUstate(dTdU, Ustate, Uwork, nspecies, ndim, mix);
        nondimensionalize_dT_dUstate(dTdU, (double*)uinf, nspecies, nd);
        //uinf14 = 0.0; uinf15 = 0.0; uinf16 = 0.0; uinf17 = 0.0; uinf18 = 0.0; uinf19 = 0.0; // LAPLACIAN
		uinf14 = dTdU[0]; uinf15 = dTdU[1]; uinf16 = dTdU[2]; uinf17 = dTdU[3]; uinf18 = dTdU[4]; uinf19 = dTdU[5];

		T t2 = udg1*udg1;
		T t3 = udg2*udg2;
		T t4 = udg3*udg3;
		T t5 = udg4*udg4;
		T t6 = udg5*udg5;
		T t7 = udg6*udg6;
		T t8 = udg7*udg7;
		T t9 = 1.0/3.141592653589793;
		T t10 = odg1*param21*uinf20;
		T t11 = 1.0/param19;
		T t12 = 1.0/param20;
		T t15 = udg1*1.0E+12;
		T t16 = udg2*1.0E+12;
		T t17 = udg3*1.0E+12;
		T t18 = udg4*1.0E+12;
		T t19 = udg5*1.0E+12;
		T t13 = t11*t12*uinf13;
		T t20 = atan(t15);
		T t21 = atan(t16);
		T t22 = atan(t17);
		T t23 = atan(t18);
		T t24 = atan(t19);
		T t40 = t2*1.0E+24;
		T t41 = t3*1.0E+24;
		T t42 = t4*1.0E+24;
		T t43 = t5*1.0E+24;
		T t44 = t6*1.0E+24;
		T t14 = t10+t13;
		T t25 = t9*t20;
		T t26 = t9*t21;
		T t27 = t9*t22;
		T t28 = t9*t23;
		T t29 = t9*t24;
		T t45 = t40+1.0;
		T t46 = t41+1.0;
		T t47 = t42+1.0;
		T t48 = t43+1.0;
		T t49 = t44+1.0;
		T t30 = t25+1.0/2.0;
		T t31 = t26+1.0/2.0;
		T t32 = t27+1.0/2.0;
		T t33 = t28+1.0/2.0;
		T t34 = t29+1.0/2.0;
		T t50 = 1.0/t45;
		T t51 = 1.0/t46;
		T t52 = 1.0/t47;
		T t53 = 1.0/t48;
		T t54 = 1.0/t49;
		T t35 = t30*udg1;
		T t36 = t31*udg2;
		T t37 = t32*udg3;
		T t38 = t33*udg4;
		T t39 = t34*udg5;
		T t60 = t9*t15*t50;
		T t61 = t9*t16*t51;
		T t62 = t9*t17*t52;
		T t63 = t9*t18*t53;
		T t64 = t9*t19*t54;
		T t55 = t35+3.182454300088011E-13;
		T t56 = t36+3.182454300088011E-13;
		T t57 = t37+3.182454300088011E-13;
		T t58 = t38+3.182454300088011E-13;
		T t59 = t39+3.182454300088011E-13;
		T t65 = t30+t60;
		T t66 = t31+t61;
		T t67 = t32+t62;
		T t68 = t33+t63;
		T t69 = t34+t64;
		T t80 = t35+t36+t37+t38+t39+1.591227150044006E-12;
		T t70 = t65*udg9;
		T t71 = t66*udg10;
		T t72 = t67*udg11;
		T t73 = t65*udg17;
		T t74 = t68*udg12;
		T t75 = t66*udg18;
		T t76 = t69*udg13;
		T t77 = t67*udg19;
		T t78 = t68*udg20;
		T t79 = t69*udg21;
		T t81 = 1.0/t80;
		T t82 = t81*t81;
		T t83 = t81*udg8;
		T t84 = t81*uinf1;
		T t85 = t81*udg6*udg7;
		T t88 = t70*t80;
		T t89 = t71*t80;
		T t90 = t72*t80;
		T t91 = t73*t80;
		T t92 = t74*t80;
		T t93 = t75*t80;
		T t94 = t76*t80;
		T t95 = t77*t80;
		T t96 = t78*t80;
		T t97 = t79*t80;
		T t110 = t70+t71+t72+t74+t76;
		T t111 = t73+t75+t77+t78+t79;
		T t86 = (t7*t82)/2.0;
		T t87 = (t8*t82)/2.0;
		T t98 = -t88;
		T t99 = -t89;
		T t100 = -t90;
		T t101 = -t91;
		T t102 = -t92;
		T t103 = -t93;
		T t104 = -t94;
		T t105 = -t95;
		T t106 = -t96;
		T t107 = -t97;
		T t108 = t83+t84;
		T t112 = t55*t110;
		T t113 = t56*t110;
		T t114 = t57*t110;
		T t115 = t58*t110;
		T t116 = t59*t110;
		T t117 = t55*t111;
		T t118 = t56*t111;
		T t119 = t57*t111;
		T t120 = t58*t111;
		T t121 = t59*t111;
		T t122 = t81*t110*udg6;
		T t123 = t81*t110*udg7;
		T t124 = t81*t111*udg6;
		T t125 = t81*t111*udg7;
		T t109 = t86+t87;
		T t126 = -t122;
		T t127 = -t123;
		T t128 = -t124;
		T t129 = -t125;
		T t134 = t98+t112;
		T t135 = t99+t113;
		T t136 = t100+t114;
		T t137 = t102+t115;
		T t138 = t104+t116;
		T t139 = t101+t117;
		T t140 = t103+t118;
		T t141 = t105+t119;
		T t142 = t106+t120;
		T t143 = t107+t121;
		T t144 = -t81*(t122-udg14);
		T t145 = -t81*(t123-udg15);
		T t146 = -t81*(t124-udg22);
		T t147 = -t81*(t125-udg23);
		T t148 = t81*(t122-udg14)*-2.0;
		T t149 = t81*(t125-udg23)*-2.0;
		T t150 = t81*(t125-udg23);
		T t152 = -t81*uinf2*(t88-t112);
		T t153 = -t82*uinf2*(t88-t112);
		T t154 = -t81*uinf3*(t89-t113);
		T t155 = -t82*uinf3*(t89-t113);
		T t156 = -t81*uinf4*(t90-t114);
		T t157 = -t82*uinf4*(t90-t114);
		T t158 = -t81*uinf5*(t92-t115);
		T t159 = -t82*uinf5*(t92-t115);
		T t160 = -t81*uinf6*(t94-t116);
		T t161 = -t82*uinf6*(t94-t116);
		T t162 = -t81*uinf2*(t91-t117);
		T t163 = -t82*uinf2*(t91-t117);
		T t164 = -t81*uinf3*(t93-t118);
		T t165 = -t82*uinf3*(t93-t118);
		T t166 = -t81*uinf4*(t95-t119);
		T t167 = -t82*uinf4*(t95-t119);
		T t168 = -t81*uinf5*(t96-t120);
		T t169 = -t82*uinf5*(t96-t120);
		T t170 = -t81*uinf6*(t97-t121);
		T t171 = -t82*uinf6*(t97-t121);
		T t172 = t81*uinf2*(t88-t112);
		T t173 = t81*uinf3*(t89-t113);
		T t174 = t81*uinf4*(t90-t114);
		T t175 = t81*uinf5*(t92-t115);
		T t176 = t81*uinf6*(t94-t116);
		T t177 = t81*uinf2*(t91-t117);
		T t178 = t81*uinf3*(t93-t118);
		T t179 = t81*uinf4*(t95-t119);
		T t180 = t81*uinf5*(t96-t120);
		T t181 = t81*uinf6*(t97-t121);
		T t187 = -t12*uinf12*(t81*(t123-udg15)+t81*(t124-udg22));
		T t194 = -t55*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116));
		T t195 = -t56*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116));
		T t196 = -t57*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116));
		T t197 = -t58*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116));
		T t198 = -t59*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116));
		T t199 = -t55*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121));
		T t200 = -t56*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121));
		T t201 = -t57*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121));
		T t202 = -t58*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121));
		T t203 = -t59*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121));
		T t130 = t126+udg14;
		T t131 = t127+udg15;
		T t132 = t128+udg22;
		T t133 = t129+udg23;
		T t151 = t150*2.0;
		T t182 = t144+t147;
		T t183 = t145+t146;
		T t184 = -odg1*(t150+t81*(t122-udg14));
		T t186 = t148+t150;
		T t188 = t12*uinf12*(t149+t81*(t122-udg14))*(-2.0/3.0);
		T t189 = t12*uinf12*(t149+t81*(t122-udg14))*(2.0/3.0);
		T t190 = t12*uinf12*(t147+t81*(t122-udg14)*2.0)*(-2.0/3.0);
		T t191 = t85+t187;
		T t192 = t153+t155+t157+t159+t161;
		T t193 = t163+t165+t167+t169+t171;
		T t185 = t144+t151;
		f[0*ng+i] = t172+t194+t55*t81*udg6;
		f[1*ng+i] = t173+t195+t56*t81*udg6;
		f[2*ng+i] = t174+t196+t57*t81*udg6;
		f[3*ng+i] = t175+t197+t58*t81*udg6;
		f[4*ng+i] = t176+t198+t59*t81*udg6;
		f[5*ng+i] = t184+t190+uinf1+t7*t81;
		f[6*ng+i] = t191;
		f[7*ng+i] = t108*udg6-uinf7*(t152+t55*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116)))-uinf8*(t154+t56*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116)))-uinf9*(t156+t57*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116)))-uinf10*(t158+t58*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116)))-uinf11*(t160+t59*(t82*uinf2*(t88-t112)+t82*uinf3*(t89-t113)+t82*uinf4*(t90-t114)+t82*uinf5*(t92-t115)+t82*uinf6*(t94-t116)))+t14*(uinf19*(udg16+t80*(t82*udg6*(t122-udg14)+t82*udg7*(t123-udg15))-t109*t110)+t70*uinf14+t71*uinf15+t72*uinf16+t74*uinf17+t76*uinf18)-t81*udg6*(odg1*(t150+t81*(t122-udg14))+t12*uinf12*(t147+t81*(t122-udg14)*2.0)*(2.0/3.0))+t81*t187*udg7;
		f[8*ng+i] = t177+t199+t55*t81*udg7;
		f[9*ng+i] = t178+t200+t56*t81*udg7;
		f[10*ng+i] = t179+t201+t57*t81*udg7;
		f[11*ng+i] = t180+t202+t58*t81*udg7;
		f[12*ng+i] = t181+t203+t59*t81*udg7;
		f[13*ng+i] = t191;
		f[14*ng+i] = t184+t189+uinf1+t8*t81;
		f[15*ng+i] = t108*udg7-uinf7*(t162+t55*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121)))-uinf8*(t164+t56*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121)))-uinf9*(t166+t57*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121)))-uinf10*(t168+t58*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121)))-uinf11*(t170+t59*(t82*uinf2*(t91-t117)+t82*uinf3*(t93-t118)+t82*uinf4*(t95-t119)+t82*uinf5*(t96-t120)+t82*uinf6*(t97-t121)))+t14*(uinf19*(udg24+t80*(t82*udg6*(t124-udg22)+t82*udg7*(t125-udg23))-t109*t111)+t73*uinf14+t75*uinf15+t77*uinf16+t78*uinf17+t79*uinf18)-t81*udg7*(t188+odg1*(t150+t81*(t122-udg14)))+t81*t187*udg6;
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);