void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = i/npe;
		dstype param1 = param[0];
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		dstype xdg2 = xdg[j+npe*1+npe*ncx*k];
		dstype t2 = xdg1*xdg1;
		dstype t3 = xdg2*xdg2;
		dstype t4 = -t2;
		dstype t5 = -t3;
		dstype t6 = t2/2.0;
		dstype t7 = t3/2.0;
		dstype t8 = -t6;
		dstype t9 = -t7;
		dstype t10 = t4+t5+1.0;
		dstype t11 = exp(t10);
		dstype t12 = t8+t9+1.0/2.0;
		dstype t13 = exp(t12);
		dstype t14 = t13*xdg1*1.591549430918953E-1;
		dstype t15 = t13*xdg2*1.591549430918953E-1;
		dstype t16 = -t15;
		dstype t17 = t14+1.0;
		f[j+npe*0+npe*nce*k] = 1.0;
		f[j+npe*1+npe*nce*k] = t16+1.0;
		f[j+npe*2+npe*nce*k] = t17;
		f[j+npe*3+npe*nce*k] = pow(t15-1.0,2.0)/2.0+t2*t11*1.266514795529222E-2+t3*t11*1.266514795529222E-2-(t11*(t2*1.266514795529222E-2+t3*1.266514795529222E-2)-1.0)/(param1-1.0)+(t17*t17)/2.0;
		f[j+npe*4+npe*nce*k] = t16;
		f[j+npe*5+npe*nce*k] = t14;
		f[j+npe*6+npe*nce*k] = 0.0;
	}
}

