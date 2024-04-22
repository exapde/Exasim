function [UDG,UH]= hdgsolve(master,mesh,pde,UDG,UH,SH)
% Solve PDEs using the HDG method and Newton iteration

nsiz = max(mesh.elcon(:));
ne = size(mesh.elcon,2);
npe = master.npe;
npf = master.npf;
nfe  = size(master.perm,2);
ncu = size(UH,1);
nc = size(UDG,2);
ncq = nc - ncu;

if isempty(SH)
  SH = zeros(npe, nc, ne);
end

if ncq>0
  [MinvC, MinvE] = qequationint(master, mesh.dgnodes);
else
  MinvC = [];
  MinvE = [];
end
[AE, FE, DUDG, DUDG_DUH] = uequationint(master,mesh,pde,UDG,UH,SH,MinvC,MinvE);

if isfield(pde, "denseblock")==0
  pde.denseblock = 0;
end
 
if pde.denseblock==0
  [K, F] = assemblelinearsystem(AE, FE, mesh.elcon);  
%   mesh.f2t = mkf2t(mesh.f, mesh.t2f);
%   F1 = assemblelRHS(FE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
%   max(abs(F(:)-F1(:))) 
%   v1 = hdgmatvec(AE, F, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
%   v2 = K*F;
%   max(abs(v1(:)-v2(:)))
%   BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));
%   w1 = applyblockjacobi(BE, v1);
%   w2 = 0*w1;
%   for i = 1:nf
%     ind = (ncu*npf*(i-1)+1):ncu*npf*i; 
%     w2(:,i) = K(ind,ind)\v2(ind);
%   end
%   max(abs(w1(:)-w2(:)))
%   pause
%   x1 = K\F;  
%   max(abs(x1(:)))
%   [x2, flags, iter, rev] = hdggmres(AE, F1, BE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]), UH(:), 1000, 1e-7, 1000);  
%   max(abs(x2(:)-x1(:)))
%   r1 = F - K*x1(:);
%   r2 = F - K*x2(:);
%   [norm(r1(:)) norm(r2(:))]
%   error("check matvec");
else
  %[K, F] = assemblelinearsystem(AE, FE, mesh.elcon, mesh.f, mesh.t2f);
  %f2f = mkf2f(mesh.f, mesh.t2f);
  
  F = assembleRHS(FE, mesh.elcon);
  if isfield(mesh, "f2t")==0
    mesh.f2t = mkf2t(mesh.f, mesh.t2f);
  end
  BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
end

it   = 0;
duh  = 1e6;
NewtonTol = 1e-8;
while duh > NewtonTol && it < 10
                
    if pde.denseblock==0
        DUH = full(reshape(K\F,ncu,nsiz));   
    else
        %DUH = gmres(K, F, f2f);
        DUH = hdggmres(AE, F, BE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]), 0*UH(:), pde.GMRESrestart, pde.linearsolvertol, pde.linearsolveriter);
        DUH = full(reshape(DUH,ncu,nsiz));    
    end

    DUHE = reshape(full(DUH(:,mesh.elcon(:))),ncu*nfe*npf,ne);    
    DUDGT = 0*DUDG;        
    for i = 1:ne
      DUDGT(:,i) = DUDG(:,i) + DUDG_DUH(:,:,i)*DUHE(:,i);
    end
    
    it = it + 1;    
    fprintf('Newton iteration :  %d\n', it);
    
    % Find stepsize for damped Newton    
    UDG0 = UDG;
    UH0  = UH;    
    duh0 = norm(F(:));   
    alfa = 1;
    while 1        
        UDG(:,1:ncu,:) = UDG0(:,1:ncu,:) + alfa*reshape(DUDGT, [npe ncu ne]);                     
        UH = UH0 + alfa*DUH;   
        
        q = qequationschur(MinvC, MinvE, UDG(:,1:ncu,:), reshape(UH(:,mesh.elcon),[ncu npf*nfe ne]), SH(:,ncu+1:nc,:), pde.fc_q);    
        UDG(:,(ncu+1):nc,:) = q;          
        
%         dataout = "/Users/ngoccuongnguyen/Exasim/build/dataout/out" + num2str(it);        
%         tmp = readbin(dataout + "newton_x.bin");        
%         max(abs(tmp(:)-DUH(:)))
%         tmp = readbin(dataout + "newton_u.bin");
%         max(abs(tmp(:)-UH(:)))
%         tmp = readbin(dataout + "newton_uh.bin");
%         max(abs(tmp(:)-UH(:)))
%         tmp = readbin(dataout + "newton_udg.bin");
%         tmp = reshape(tmp, size(UDG));
%         tmu = tmp(:,1:ncu,:) - UDG(:,1:ncu,:);
%         max(abs(tmu(:)))
%         tmq = tmp(:,ncu+1:end,:) - q;
%         max(abs(tmq(:)))
%         
%         tmp = readbin(dataout + "newton_FE.bin");                
%         max(abs(tmp(:)-FE(:)))
%         
%         tmp = readbin(dataout + "newton_AE.bin");                
%         max(abs(tmp(:)-AE(:)))
%         
%         tmp = readbin(dataout + "newton_DUDG.bin");                
%         max(abs(tmp(:)-DUDG(:)))
%         
%         tmp = readbin(dataout + "newton_DUDG_DUH.bin");                
%         max(abs(tmp(:)-DUDG_DUH(:)))
%         pause
%         

        [AE, FE, DUDG, DUDG_DUH] = uequationint(master,mesh,pde,UDG,UH,SH,MinvC,MinvE);
        if pde.denseblock==0
          [K, F] = assemblelinearsystem(AE, FE, mesh.elcon);
        else
          F = assembleRHS(FE, mesh.elcon);
          BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
        end

        duh  = norm(F(:));                         
        if (duh>duh0) 
            alfa=alfa/2;
            fprintf(' alpha = %f\n',alfa);
            if alfa<1e-2, break; end
        else
            break;
        end
    end
               
    fprintf('Old residual: %e,   New residual: %e    %e\n', [duh0 duh alfa]);       
end

if ncq>0
  q = qequationschur(MinvC, MinvE, UDG(:,1:ncu,:), reshape(UH(:,mesh.elcon),[ncu npf*nfe ne]), SH(:,ncu+1:nc,:), pde.fc_q);    
  UDG(:,(ncu+1):end,:) = q;  
end

end


function f2t = mkf2t(f, t2f)

nf = size(f,1);
f2t = zeros(4,nf);
for i = 1:nf
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
        i2 = find(kf(2,:)==i);  % obtain the index of face i in the 2nd element                                                    
        f2t(1,i) = fi(1);
        f2t(2,i) = i1;
        f2t(3,i) = fi(2);
        f2t(4,i) = i2;
    else % face i is a boundary face
        kf = t2f(fi(1),:); % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element                 
        f2t(1,i) = fi(1);
        f2t(2,i) = i1;        
    end        
end 

end

