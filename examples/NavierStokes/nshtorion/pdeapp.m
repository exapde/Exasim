pdeapp_ns;
pdeapp_ht;

% axis symmetric, supersonic outflow, supersonic inflow, nonhomogenous thermal wall
mesh.boundarycondition = [6, 2, 1, 7]; % Set boundary condition for each boundary
mesh.ibwall = 4;   % interface boundary index for fluid domain

% axis symmetric (homogenuous Neumann), heat flux (nonhomogenuous Neumann), iso-thermal (nonhomogenuous Dirichlet) 
meshht.boundarycondition = [3, 4, 1]; % Set boundary condition for each boundary
meshht.ibwall = 2; % interface boundary index for solid domain

iter = 0; itermax = 10;
while (iter < itermax)
  iter = iter + 1;
    
  if iter>1
    meshht.udg = solht;
  end
  % transfer heat flux from fluid domain to solid domain
  [meshht.vdg, UHbns] = solutiontransfer_ns2ht(pde, dmd, mesh, meshht);
  
  % solve the heat equation with the external heat flux from fluid domain
  if iter == 1
      [solht,pdeht,meshht,masterht,dmdht] = exasim(pdeht,meshht);
  else
    [pdeht,meshht,masterht,dmdht] = preprocessing(pdeht,meshht);
    runcode(pdeht, 1); 
    solht = fetchsolution(pdeht,masterht,dmdht, pdeht.buildpath + '/dataout');
  end 
  figure(1); clf; scaplot(meshht, solht(:,1,:));
  set(gca,'FontSize',20); axis equal; axis tight;

  % transfer temperature from solid domain to fluid domain
  mesh.udg = sol;
  [mesh.vdg(:,2,:), UHbht, in, im] = solutiontransfer_ht2ns(pdeht, dmdht, meshht, mesh);

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
  
  % solve the NS equations with the external temperature from solid domain
  [pde,mesh,master,dmd] = preprocessing(pde,mesh);
  runcode(pde, 1); % run C++ code
  sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
  figure(2); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
  set(gca,'FontSize',20); axis equal; axis tight;
  figure(3); clf; scaplot(mesh, eulereval(sol, 't',gam,Minf),[]);
  set(gca,'FontSize',20); axis equal; axis tight;
end

solht = fetchsolution(pdeht,masterht,dmdht, pdeht.buildpath + '/dataout');
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(4); clf; scaplot(meshht, (Tref/Tinf)*solht(:,1,:), [200 5000]);
hold on;
scaplot(mesh, (Tref/Tinf)*eulereval(sol, 't',gam,Minf),[200 5000]);
set(gca,'FontSize',20); 
%set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
xlabel("$z$", 'interpreter', 'latex', 'FontSize', 28);
ylabel("$r$", 'interpreter', 'latex', 'FontSize', 28);



