mesh.porder = pde.porder;
it = 2400
tmp = getsolution(['dataout/out_t' num2str(it)],dmd,master.npe);
p = eulereval(tmp, "p", gam, Minf);
figure(1); clf; scaplot(mesh,p,[])
