function [AE, FE, DUDG, DUDG_DUH] = uequationint(master,mesh,pde,UDG,UH,SH,MinvC,MinvE)

ne  = size(mesh.dgnodes,3);
ns  = 128;
nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    

nd = master.nd;
pde.nc = size(UDG,2);
pde.ncu = size(UH,1);
pde.ncq = pde.nc - pde.ncu;
npe = master.npe;
npf = master.npf;
nfe = size(master.perm,2);
ndf = size(mesh.elcon,1);
ncu = pde.ncu;
nc = pde.nc;
ncq = nc - ncu;

if isfield(pde, 'debugmode') == 0  
  pde.debugmode = 0;
end
if pde.debugmode == 1
  if ncq>0
    save qequationint.mat MinvC MinvE;
  end
end
    
if isempty(SH)
  SH = zeros(npe, nc, ne);
end

UHAT = UH;
UH  = reshape(UH(:,mesh.elcon),[ncu ndf ne]);

% solve q + nabla u = 0 or dq/dt + nabla u = 0 to obtain q
if (ncq > 0)  
  q = qequationschur(MinvC, MinvE, UDG(:,1:ncu,:), UH, SH(:,ncu+1:nc,:), pde.fc_q);    
  UDG(:,(ncu+1):nc,:) = q;  
  % MinvC = npe * npe * ne * nd
  MinvC = reshape(permute(MinvC, [1 4 2 3]),[npe*nd npe ne]);
  % MinvE = npe * ndf * ne * nd
  MinvE = reshape(permute(MinvE, [1 4 2 3]),[npe*nd ndf ne]);        
end

if pde.debugmode == 1      
  save uequationsol.mat UDG UHAT;
end

UH = permute(UH,[2 3 1]);
UDG = permute(UDG,[1 3 2]);
SH  = permute(SH,[1 3 2]);
dgnodes = permute(mesh.dgnodes,[1 3 2]);

DUDG = zeros(npe*ncu, ne);
DUDG_DUH = zeros(npe*ncu, ncu*ndf, ne);
FE = zeros(ncu*ndf, ne);
AE = zeros(ncu*ndf, ncu*ndf, ne);
for j=1:nb  
    id = nm(1,j):nm(2,j);        
        
    % compute (du/dt, w), (Flux, nabla w) and (Source, w)
    [Ru, BD, pg, udgg, f, f_udg, s, s_udg] = uequationelemint(master,pde,dgnodes(:,id,:),UDG(:,id,:),SH(:,id,:));            
    
    if pde.debugmode == 1      
      save uequationelemint.mat Ru BD pg udgg f f_udg s s_udg;
    end
    
    % compute <Fhat, w> and <Fbou, mu> 
    [Ru, Rh, B, D, F, G, K, H] = uequationfaceint(master,pde,mesh.bf(:,id),dgnodes(:,id,:),UDG(:,id,:),UH(:,id,:),Ru,BD);
    if pde.debugmode == 1
      save uequationfaceint.mat Ru Rh B D F G K H;
    end
    [norm(Ru(:)) norm(B(:)) norm(D(:)) norm(Rh(:))] 
    [norm(F(:)) norm(K(:)) norm(G(:)) norm(H(:))]
    
    % solve the local problem 
    if ncq > 0
      [AE(:,:,id), FE(:,id), DUDG(:,id), DUDG_DUH(:,:,id), D, F, K, H] = uequationschur(Ru, Rh, B, D, F, G, K, H, MinvC(:,:,id), MinvE(:,:,id));
    else
      [AE(:,:,id), FE(:,id), DUDG(:,id), DUDG_DUH(:,:,id), D, F, K, H] = uequationschur(Ru, Rh, B, D, F, G, K, H);
    end
    t1 = DUDG(:,id);
    t2 = FE(:,id);
    t3 = DUDG_DUH(:,:,id);
    t4 = AE(:,:,id);
    [norm(D(:)) norm(F(:)) norm(K(:)) norm(H(:))]
    [norm(t1(:)) norm(t2(:)) norm(t3(:)) norm(t4(:))]
    
    if pde.debugmode == 1
      save uequationschur.mat AE FE DUDG DUDG_DUH D F K H;            
      R = assemblelRHS(FE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
      save assemblelRHS.mat R      
      v = hdgmatvec(AE, R, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
      save matvec.mat v           
      BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
      save blockjacobi.mat BE        
      w = applyblockjacobi(BE, v);
      save applyblockjacobi.mat w  
      error("Stop for debugging Exasim C++ ...");
    end    
end    


