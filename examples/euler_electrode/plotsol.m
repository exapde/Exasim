udg=getsolution("/Users/cuongnguyen/Documents/GitHub/Exasim/build/dataout/out",dmd,6);
r=udg(:,1,:);
rmin=min(r(:));
i=find(r(:)==rmin);
x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
figure(1); clf; scaplot(mesh,udg(:,1,:),[],1); hold on; plot(x(i),y(i),'or');axis on; axis equal; axis tight;

n = 50;
[k,sol] = readsoldmd("/Users/cuongnguyen/Documents/GitHub/Exasim/build/dataout/outudg", dmd, n, 0);
figure(2); clf; scaplot(mesh,sol(:,2,:,end)./sol(:,1,:,end),[],1); hold on; plot(x(i),y(i),'or');axis on; axis equal; axis tight;
figure(3); clf; scaplot(mesh,sol(:,3,:,end)./sol(:,1,:,end),[],1); hold on; plot(x(i),y(i),'or');axis on; axis equal; axis tight;
figure(4); clf; scaplot(mesh,sol(:,1,:,end),[],1); hold on; plot(x(i),y(i),'or');axis on; axis equal; axis tight;
