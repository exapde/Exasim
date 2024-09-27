dRu_AD = getsolution('dataout/outenz_dRu',dmd,master.npe);
dRu_FD = getsolution('dataout/outenz_fd',dmd,master.npe);
% u_AD = getsolution('dataout/outenz_udg',dmd,master.npe);
figure(1); clf; hold on; plot(dRu_AD(:),'-o'); plot(dRu_FD(:))
figure(2); clf; plot(dRu_FD(:));
diff = (dRu_AD(:) - dRu_FD(:)); %./ (dRu_FD(:) + 1e-8);
figure(3); clf; plot(diff)
