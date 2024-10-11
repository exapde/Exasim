fname = "/home/austinsp/Exasim/build/dataout/out_t1100";
UDG = getsolution(fname, dmd, 6);

normE = sqrt(UDG(:,6,:).^2+UDG(:,9,:).^2);

scaplot(mesh, normE,[],0,1); colormap jet;