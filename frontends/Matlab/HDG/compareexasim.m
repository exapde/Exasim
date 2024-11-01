function compareexasim(master, mesh, pde)

e = zeros(32,1);

load qequationint.mat;
load uequationsol.mat;
load uequationelemint.mat;

ncu = size(Ru,3);
nc = size(BD,5);
pde.ncu = ncu;
pde.nc = nc;
pde.ncq = nc - ncu;
nge = master.nge;
npe = master.npe;
npf = master.npf;
nd = master.nd;
nfe = size(master.perm,2);
ne = size(mesh.t, 2);
nf = size(mesh.f2t,2);

D = BD(:,:,:,:,1:ncu);
B = BD(:,:,:,:,(ncu+1):end);

dataout = pde.buildpath + "/dataout/out";

[MinvC1, MinvE1, ~, ~, uqg1, uhat1] = qEquationExasim(dataout, pde, npe, npf, nd, ne);
[Ru1, B1, D1, fg1, fg_udg1, sg1, sg_udg1, xg1, udg1] = uEquationElemExasim(dataout, pde, nge, npe, ne);

e(1) = max(abs(MinvC(:)-MinvC1(:)));
e(2) = max(abs(MinvE(:)-MinvE1(:)));
e(3) = max(abs(Ru(:)-Ru1(:)));
e(4) = max(abs(B(:)-B1(:)));
e(5) = max(abs(D(:)-D1(:)));

[Ru1, Rh1, B1, D1, F1, G1, K1, H1] = uEquationFaceExasim(dataout, pde, npe, npf, nfe, ne);
load uequationfaceint.mat;

e(6) = max(abs(Rh(:)-Rh1(:)));
e(7) = max(abs(Ru(:)-Ru1(:)));
e(8) = max(abs(B(:)-B1(:)));
e(9) = max(abs(D(:)-D1(:)));
e(10) = max(abs(F(:)-F1(:)));
e(11) = max(abs(G(:)-G1(:)));
e(12) = max(abs(K(:)-K1(:)));
e(13) = max(abs(H(:)-H1(:)));

[AE1, FE1, DUDG1, DUDG_DUH1, D1, F1, K1, H1] = uEquationSchurExasim(dataout, pde, npe, npf, nfe, ne);
load uequationschur.mat;

e(14) = max(abs(D(:)-D1(:)));
e(15) = max(abs(F(:)-F1(:)));
e(16) = max(abs(K(:)-K1(:)));
e(17) = max(abs(H(:)-H1(:)));
e(18) = max(abs(DUDG(:)-DUDG1(:)));
e(19) = max(abs(DUDG_DUH(:)-DUDG_DUH1(:)));
e(20) = max(abs(FE(:)-FE1(:)));
e(21) = max(abs(AE(:)-AE1(:)));

R1 = assembleRhsExasim(dataout, pde, npf, nf);
load assemblelRHS.mat;
e(22) = max(abs(R(:)-R1(:)));

BE1 = blockjacobiExasim(dataout, pde, npf, nf);
load blockjacobi.mat;
e(23) = max(abs(BE(:)-BE1(:)));

% v1 = matvecExasim(pde, mesh);
% load matvec.mat;
% e(24) = max(abs(v(:)-v1(:)));

disp("Check udg and xdg...");
pg = pg(:,1:nd); % For matlab, AV continuation uses dgnodes to store av field
eu = UDG(:,1:ncu,:)-uqg1(:,1:ncu,:);
fprintf('Maximum Absolue Error in U =  %e \n', max(abs(eu(:))));       
if (nc>ncu)
  eq = abs(UDG(:,ncu+1:nc,:)-uqg1(:,ncu+1:nc,:));
  fprintf('Maximum Absolue Error in Q =  %e \n', max(abs(eq(:))));       
end
fprintf('Maximum Absolue Error in UHAT =  %e \n', max(abs(UHAT(:)-uhat1(:))));    
fprintf('Maximum Absolue Error in UDG at Gauss points =  %e \n', max(abs(udgg(:)-udg1(:))));       
fprintf('Maximum Absolue Error in dgnodes =  %e \n', max(abs(pg(:)-xg1(:))));       
fprintf('\n');

disp("Check qEquation...");
fprintf('Maximum Absolue Error in MinvC =  %e \n', e(1));       
fprintf('Maximum Absolue Error in MinvE =  %e \n', e(2));       
fprintf('\n');

disp("Check source term and flux...");
f_udg = permute(f_udg, [1 2 3 5 4]);
fprintf('Maximum Absolue Error in source term =  %e \n', max(abs(s(:)-sg1(:))));   
fprintf('Maximum Absolue Error in source term derivative =  %e \n', max(abs(s_udg(:)-sg_udg1(:))));   
fprintf('Maximum Absolue Error in flux =  %e \n', max(abs(f(:)-fg1(:))));   
fprintf('Maximum Absolue Error in flux derivative =  %e \n', max(abs(f_udg(:)-fg_udg1(:))));   
fprintf('\n');

disp("Check uEquationElem...");
fprintf('Maximum Absolue Error in Ru =  %e \n', e(3));       
fprintf('Maximum Absolue Error in B =  %e \n', e(4));       
fprintf('Maximum Absolue Error in D =  %e \n', e(5));   
fprintf('\n');

disp("Check uEquationFace...");
fprintf('Maximum Absolue Error in Rh =  %e \n', e(6));      
fprintf('Maximum Absolue Error in Ru =  %e \n', e(7));       
fprintf('Maximum Absolue Error in B =  %e \n', e(8));       
fprintf('Maximum Absolue Error in D =  %e \n', e(9));    
fprintf('Maximum Absolue Error in F =  %e \n', e(10));    
fprintf('Maximum Absolue Error in G =  %e \n', e(11));    
fprintf('Maximum Absolue Error in K =  %e \n', e(12));    
fprintf('Maximum Absolue Error in H =  %e \n', e(13));    
fprintf('\n');

disp("Check uEquationSchur...");
fprintf('Maximum Absolue Error in D =  %e \n', e(14));    
fprintf('Maximum Absolue Error in F =  %e \n', e(15));    
fprintf('Maximum Absolue Error in K =  %e \n', e(16));    
fprintf('Maximum Absolue Error in H =  %e \n', e(17));    
fprintf('Maximum Absolue Error in DUDG =  %e \n', e(18));    
fprintf('Maximum Absolue Error in DUDG_DUH =  %e \n', e(19));    
fprintf('Maximum Absolue Error in FE =  %e \n', e(20));    
fprintf('Maximum Absolue Error in AE =  %e \n', e(21));    

disp("Check assembleRHS...");
fprintf('Maximum Absolue Error in RHS =  %e \n', e(22));    

disp("Check blockjacobi...");
fprintf('Maximum Absolue Error in block jacobi =  %e \n', e(23));    

% disp("Check matvec...");
% fprintf('Maximum Absolue Error in MatVec =  %e \n', e(24));    


