% pdeapp_ns; 
% pdeapp_ht;

%% Uncomment these lines if you want to use the r-adaptive solution
pdeapp_ns_adapt; mesh1 = mesh_adapt;
pdeapp_ht_adapt; meshht = meshht2;

mesh1.boundarycondition   = [11;5;9]; 
pde_ns.bcm                = [11,5,9];  
mesh1.ibwall = 1;

meshht.boundarycondition = [4;3;1]; 
meshht.ibwall = 1; 
%%
iter = 0; itermax = 10;
while (iter < itermax)
  iter = iter + 1;
  
  % transfer heat flux from fluid domain to solid domain
  meshht.udg = solht;
  [meshht.vdg, UHbns] = solutiontransfer_ns2ht(pde_ns, dmd, mesh1, meshht, UDG0, UH0);
  
  % solve the heat equation with the external heat flux from fluid domain
  [pdeht,meshht,masterht,dmdht] = preprocessing(pdeht,meshht);
  runcode(pdeht, 1); 
  solht = fetchsolution(pdeht,masterht,dmdht, pdeht.buildpath + '/dataout');
  figure(1); clf; scaplot(meshht, solht(:,1,:));
  set(gca,'FontSize',20); axis equal; axis tight;

  % transfer temperature from solid domain to fluid domain
  % mesh.udg = sol;
  
  [tmp_ht2ns, UHbht, in, im] = solutiontransfer_ht2ns(pdeht, dmdht, meshht, mesh1);
  mesh1.dgnodes(:,4,:) = tmp_ht2ns;

  % check convergence if the solid temperature matches the fluid temperature
  Tns = eulereval(UHbns, 't', gam, Minf);
  Tht = UHbht;
  dT = 0*Tht;
  for i = 1:length(in)
    i2 = in(i);
    dT(:,:,i) = Tht(:,:,i) - Tns(im(:,i),:,i2) ;
  end  
  fprintf('|T_f - T_s| = %f\n', norm(dT(:)));
  if norm(dT(:)) < 1e-3
    break;
  end
  pause(1)

  [UDG,UH] = hdgsolve(master,mesh1,pde_ns,UDG0,UH0,[]);
  UDG0 = UDG; UH0 = UH;
end
%%
mesh1.dgnodes = mesh1.dgnodes(:,1:2,:);
mesh1.xpe = mesh1.plocal;
mesh1.telem= master.telem;
solht = fetchsolution(pdeht,masterht,dmdht, pdeht.buildpath + '/dataout');
figure(4); clf; scaplot(meshht, (Tref/Tinf)*solht(:,1,:));
hold on;
scaplot(mesh1, (Tref/Tinf)*eulereval(UDG, 't',gam,Minf),[]);
set(gca,'FontSize',20); axis equal; axis tight;


