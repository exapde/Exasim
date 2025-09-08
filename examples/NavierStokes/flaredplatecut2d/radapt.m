
av = mesh.vdg(:,1,:);
av = av/max(av(:));
rho = 1 + 3*av;

grad = sqrt(sol(:,5,:).^2 + sol(:,9,:).^2)./sol(:,1,:);
grad = divergence(sol, 1);
grad = grad/max(grad(:));
grad = limiting(grad, 0, 0.25, 1e2, 0);
dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,[2]); % distance to the wall
rho = 1 + 40*grad.*tanh(1e2*dist);
figure(2); clf; scaplot(mesh, rho,[],2);

solhm = hmsmoothing(mesh, rho, 0.1, 10);
mesh1 = radaptivity(mesh, solhm(:,1,:), -solhm(:,2,:), -solhm(:,3,:), 1e-6);


