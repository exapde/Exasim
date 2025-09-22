function [UDG,UH,AE,FE,K,F]= hdgsolve(master,mesh,pde,UDG,UH,SH)
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
  if pde.denseblock==1
    BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
  elseif pde.denseblock==2
    BE = elementaladditiveschwarz(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
  elseif pde.denseblock==3   
    [mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf] = pathreordering(mesh.epath, mesh.nelems, mesh.t2f);
    [BE.A, BE.B, BE.C, BE.D] = pathcompute(AE, mesh.epath, mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf, mesh.nelems, mesh.f2t, mesh.t2f, reshape(mesh.elcon, [npf nfe ne]));
    BE.nelems = mesh.nelems;
    BE.fpath = mesh.fpath;
    BE.fintf = mesh.fintf;
  elseif pde.denseblock==4   
    [mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf] = pathreordering2(mesh.epath, mesh.t2f);
    [BE.A, BE.B1, BE.B2, BE.C1, BE.C2, BE.D1, BE.D2, BE.DL, BE.DU] = pathcompute2(AE, mesh.epath, mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf, mesh.f2t, mesh.t2f, reshape(mesh.elcon, [npf nfe ne]));
    BE.fpath = mesh.fpath;
    BE.fintf = mesh.fintf;
    BE.nep = size(mesh.epath,2);
  elseif pde.denseblock==5
    [mesh.row_ptr, mesh.col_ind, mesh.face] = facereordering(mesh.elem, mesh.t, mesh.t2f, mesh.elemtype, size(mesh.p,1));    
    %[mesh.row_ptr, mesh.col_ind, mesh.face] = crs_faceordering(mesh.elem, mesh.f2t);    
    [neb, nfeb] = size(mesh.face);
    BE.A = zeros(ncu*npf*nfeb, ncu*npf*nfeb, neb);
    for i = 1:neb
      [~, BE.A(:,:,i)] = crs_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, mesh.face(i,:), mesh.row_ptr, mesh.col_ind, ncu, npf, nfe);
    end
    BE.face = mesh.face;
  elseif pde.denseblock==6
    %[BE.row_ptr, BE.col_ind, BE.face] = facereordering(mesh.elem, mesh.t, mesh.t2f, mesh.elemtype, size(mesh.p,1));    
    [BE.row_ptr, BE.col_ind, BE.face] = crs_faceordering(mesh.elem, mesh.f2t);    
    [ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl, BE.Lind_ji, BE.Uind_ji, BE.Lnum_ji, BE.Unum_ji] = crs_indexingilu0(BE.row_ptr, BE.col_ind, nfe);
%     max(abs(BE.face(:)-face(:)))
%     max(abs(BE.row_ptr(:)-row_ptr(:)))
%     max(abs(BE.col_ind(:)-col_ind(:)))
%     BE.A = zeros(ncu*npf, ncu*npf, BE.row_ptr(end), neb);
%     for i = 1:neb
%       val = crs_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face(i,:), BE.row_ptr, BE.col_ind, ncu, npf, nfe);      
%       BE.A(:,:,:,i) = block_ilu0(BE.row_ptr, BE.col_ind, val);       
%     end
    BE.A = crs_fullassembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face, BE.row_ptr, BE.col_ind, ncu, npf, nfe);      
    %BE.A = crs_blockilu0(BE.row_ptr, BE.col_ind, BE.A);       
    BE.A = crs_parblockilu0(ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl, BE.A);
  elseif pde.denseblock==7
    [BE.row_ptr, BE.col_ind, BE.face, ~, BE.row_ptr2, BE.col_ind2, ~, BE.color, BE.idr1, BE.idr2, BE.idr3, BE.idx1, BE.idx2, BE.idx3] = facereordering(mesh.elem, mesh.t, mesh.t2f, mesh.elemtype, size(mesh.p,1));    
    count1 = size(BE.idr1,2);
    count2 = size(BE.idr2,2);
    count3 = size(BE.idr3,2);
    count4 = BE.row_ptr2(end);
    [BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, BE.C1, BE.C2, BE.C3, BE.D] = matrix_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face, BE.row_ptr, BE.col_ind, BE.color, count1, count2, count3, count4, ncu, npf, nfe);
    [BE.A1, BE.A2, BE.A3, BE.C1, BE.C2, BE.C3, BE.D] = matrix_compute(BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, BE.C1, BE.C2, BE.C3, BE.D, BE.idx1, BE.idx2, BE.idx3);
    [neb, ~] = size(BE.face);
    for i = 1:neb
       BD = block_ilu0(BE.row_ptr2, BE.col_ind2, squeeze(BE.D(:,:,i,:)));   
       BE.D(:,:,i,:) = reshape(BD, [npf*ncu npf*ncu 1 size(BD,3)]);
    end
  end
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
          if pde.denseblock==1
            BE = blockjacobi(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
          elseif pde.denseblock==2
            BE = elementaladditiveschwarz(AE, mesh.f2t, reshape(mesh.elcon, [npf nfe ne]));  
          elseif pde.denseblock==3          
            [BE.A, BE.B, BE.C, BE.D] = pathcompute(AE, mesh.epath, mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf, mesh.nelems, mesh.f2t, mesh.t2f, reshape(mesh.elcon, [npf nfe ne]));  
          elseif pde.denseblock==4          
            [BE.A, BE.B1, BE.B2, BE.C1, BE.C2, BE.D1, BE.D2, BE.DL, BE.DU] = pathcompute2(AE, mesh.epath, mesh.fpath, mesh.lpath, mesh.fintf, mesh.lintf, mesh.f2t, mesh.t2f, reshape(mesh.elcon, [npf nfe ne]));
          elseif pde.denseblock==5
            for i = 1:neb
              [~, BE.A(:,:,i)] = crs_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, mesh.face(i,:), mesh.row_ptr, mesh.col_ind, ncu, npf, nfe);
            end  
          elseif pde.denseblock==6
%             for i = 1:neb
%               val = crs_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face(i,:), BE.row_ptr, BE.col_ind, ncu, npf, nfe);      
%               BE.A(:,:,:,i) = block_ilu0(BE.row_ptr, BE.col_ind, val);       
%             end
            BE.A = crs_fullassembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face, BE.row_ptr, BE.col_ind, ncu, npf, nfe);      
            %BE.A = crs_blockilu0(BE.row_ptr, BE.col_ind, BE.A);       
            BE.A = crs_parblockilu0(ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl, BE.A);
          elseif pde.denseblock==7
            [BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, BE.C1, BE.C2, BE.C3, BE.D] = matrix_assembly(AE, reshape(mesh.elcon, [npf nfe ne]), mesh.f2t, BE.face, BE.row_ptr, BE.col_ind, BE.color, count1, count2, count3, count4, ncu, npf, nfe);
            [BE.A1, BE.A2, BE.A3, BE.C1, BE.C2, BE.C3, BE.D] = matrix_compute(BE.A1, BE.A2, BE.A3, BE.B1, BE.B2, BE.B3, BE.C1, BE.C2, BE.C3, BE.D, BE.idx1, BE.idx2, BE.idx3);
            [neb, ~] = size(BE.face);
            for i = 1:neb
               BD = block_ilu0(BE.row_ptr2, BE.col_ind2, squeeze(BE.D(:,:,i,:)));   
               BE.D(:,:,i,:) = reshape(BD, [npf*ncu npf*ncu 1 size(BD,3)]);
            end          
          end
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

