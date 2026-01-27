% runstr = '!python3 move_files_from_hpc.py';
% eval(char(runstr));

mesh1 = hdgmesh(mesh,pde.porder);


param = {};
for ii = 1:12
    param{ii} = pde.externalparam(ii);
end
pde.arg = param;

for ii = iplt
    UDG_plt = UDG1{ii};
    UH_plt = UH1{ii};
    mesh1.dgnodes(:,3,:) =0;

[UDG_dim, UH_dim] = UDG_nondim_to_dim(param, UDG_plt, UH_plt);

[Cp,~,x1,~,~,~,Ch,~] = getsurfacedata(master, mesh1, pde, UDG_dim, UH_dim, 1);
theta1 = atan2(x1(:,2),-x1(:,1))*180/pi;
C = colororder;
figure(100); hold on; plot(theta1, -Cp,drawstr,"LineWidth",2)
figure(105); hold on; plot(theta1, Ch,drawstr,"LineWidth",2); 

end