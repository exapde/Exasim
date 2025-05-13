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
		dstype t6 = t2*(1.0/2.0);
		dstype t7 = t3*(1.0/2.0);
		dstype t4 = -t6-t7+1.0/2.0;
		dstype t5 = exp(t4);
		dstype t8 = t5*xdg1*1.591549430918953E-1;
		dstype t9 = t8+1.0;
		dstype t13 = t5*xdg2*1.591549430918953E-1;
		dstype t10 = t13-1.0;
		dstype t11 = -t2-t3+1.0;
		dstype t12 = exp(t11);
		f[j+npe*0+npe*nce*k] = 1.0;
		f[j+npe*1+npe*nce*k] = t5*xdg2*(-1.591549430918953E-1)+1.0;
		f[j+npe*2+npe*nce*k] = t9;
		f[j+npe*3+npe*nce*k] = t2*t12*1.266514795529222E-2+t3*t12*1.266514795529222E-2-(t12*(t2*1.266514795529222E-2+t3*1.266514795529222E-2)-1.0)/(param1-1.0)+(t9*t9)*(1.0/2.0)+(t10*t10)*(1.0/2.0);
		f[j+npe*4+npe*nce*k] = -t13;
		f[j+npe*5+npe*nce*k] = t8;
		f[j+npe*6+npe*nce*k] = 0.0;
	}
}

