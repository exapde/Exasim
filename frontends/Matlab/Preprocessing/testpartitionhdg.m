mpiprocs = 4;
n = 4;
elemtype = 1;
nd = 2;
porder = 2;

[mesh.p,mesh.t] = squaremesh(4,4,1,elemtype);
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;2;3;4]; 

[mesh.f, mesh.tprd, t2t, fprd] = facenumbering(mesh.p,mesh.t,elemtype,mesh.boundaryexpr,mesh.periodicexpr);
f2t = mkf2e(mesh.t,elemtype,nd);
t2f = mke2f(f2t);
dmd = meshpartitionhdg(mesh.tprd,mesh.f,t2t,mesh.boundarycondition,nd,elemtype,porder,mpiprocs,"mpmetis");

plotpartition(mesh.p',mesh.t',mesh.f',dmd);

return;

% fix f2t to ensure DOF consistency across interface faces
% the global ID of the FIRST face is smaller than that of the SECOND face    
for i=1:length(dmd)
  N = size(dmd{i}.f2t,2);
  for face = 1:N % loop over each face
    e1 = dmd{i}.f2t(1,face);
    e2 = dmd{i}.f2t(3,face);
    if e2>0
      g1 = dmd{i}.elempart(e1);
      g2 = dmd{i}.elempart(e2);
      if (g2<g1)      
        tm = dmd{i}.f2t(3:4,face);
        dmd{i}.f2t(3:4,face) = dmd{i}.f2t(1:2,face);
        dmd{i}.f2t(1:2,face) = tm;            
      end
      e1 = dmd{i}.f2t(1,face);
      e2 = dmd{i}.f2t(3,face);
      g1 = dmd{i}.elempart(e1);
      g2 = dmd{i}.elempart(e2);
      if (g2<g1)
        error("something wrong");
      end
    end
  end  
end


% fix f2t to ensure DOF consistency across interface faces
% the global ID of the FIRST face is smaller than that of the SECOND face    
for i=1:length(dmd)
  N = length(dmd{i}.elemsend);
  for n = 1:N % loop over each interface element/face
    esend = dmd{i}.elemsend(n);
    erecv = dmd{i}.elemrecv(n);
    for j = 1:nfe
      fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
      for k = 1:nfe
        frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
        if (fsend == frecv)          % if esend and erecv share the same face
          dmd{i}.facesend(n) = j;    % local ID of the face sent to neighbor
          dmd{i}.facerecv(n) = k;    % local ID of the face received from neighbor                      
          face = fsend;
        end
      end
    end                    
    e1 = dmd{i}.f2t(1,face);
    e2 = dmd{i}.f2t(3,face);
    g1 = dmd{i}.elempart(e1);
    g2 = dmd{i}.elempart(e2);
    if (g2<g1)      
      tm = dmd{i}.f2t(3:4,face);
      dmd{i}.f2t(3:4,face) = dmd{i}.f2t(1:2,face);
      dmd{i}.f2t(1:2,face) = tm;            
    end
    e1 = dmd{i}.f2t(1,face);
    e2 = dmd{i}.f2t(3,face);
    g1 = dmd{i}.elempart(e1);
    g2 = dmd{i}.elempart(e2);
    if (g2<g1)
      error("something wrong");
    end
  end  
end

% for i=1:length(dmd)
%   gsend = dmd{i}.elempart(dmd{i}.elemsend);
%   grecv = dmd{i}.elempart(dmd{i}.elemrecv);    
%   dmd{i}.gsendrecv = [gsend(:) grecv(:) 0*grecv(:)];
% end
% % 1 -> keep, 2 -> remove
% dmd{1}.gsendrecv(1:2:end,3) = 1;
% dmd{1}.gsendrecv(2:2:end,3) = 2;
% for i=2:length(dmd)
%   N = length(dmd{i}.elemsend);
%   token = 1;
%   for j = 1:N
%     nb = dmd{i}.nbsend(j);
%     if nb < i      
%       for k = 1:length(dmd{nb}.elemsend)
%         if (dmd{i}.gsendrecv(j,1) == dmd{nb}.gsendrecv(k,2)) && (dmd{i}.gsendrecv(j,2) == dmd{nb}.gsendrecv(k,1))
%           dmd{i}.gsendrecv(j,3) = 3 - dmd{nb}.gsendrecv(k,3);        
%         end
%       end
%     else
%       dmd{i}.gsendrecv(j,3) = token;
%       token = 3-token;
%     end    
%   end
% end

for i=1:length(dmd)
  N = length(dmd{i}.elemsend);
  dmd{i}.removed = zeros(N,1);  
  for n = 1:N
    if (dmd{i}.gsendrecv(n,3)==2) % face to be removed
      esend = dmd{i}.elemsend(n); % element sent to neighbor
      erecv = dmd{i}.elemrecv(n); % element received from neighbor
      for j = 1:nfe
        fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
        for k = 1:nfe
          frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
          if (fsend == frecv)          % if esend and erecv share the same face          
            dmd{i}.removed(n) = fsend;
          end
        end
      end
    end
  end
   dmd{i}.f2t(:,dmd{i}.removed(dmd{i}.removed > 0)) = [];
end

npf = porder + 1;
b = cell(4,1);
nrmb = 0;
nrmc = 0;
for i = 1:4
  dataout = pde.buildpath + "/dataout/out" + num2str(i-1);
  filename = dataout + "hdgAssembleRHS.bin";
  tmp = readbin(filename);
  b{i} = reshape(tmp, npf, length(tmp)/npf);
  nrmb = nrmb + norm(b{i}(:))^2;
  ne1 = sum(dmd{i}.elempartpts(1:2));
  for j = 1:size(b{i},2)
    bj = b{i}(:,j);
    fj = dmd{i}.f2t(:,j);
    e1 = fj(1);  
    e2 = fj(3);
    if (e1 > ne1) || (e2 > ne1)
      nrmc = nrmc + norm(bj).^2;
    end
  end
end
nrmr = sqrt(nrmb - 0.5*nrmc);
nrmb = sqrt(nrmb);
nrmc = sqrt(nrmc);
nrmR = sqrt(sum(sum(R.*R)));

npf = porder + 1;
b = cell(4,1);
for i = 1:4
  dataout = pde.buildpath + "/dataout/out" + num2str(i-1);
  filename = dataout + "hdgAssembleRHS.bin";
  tmp = readbin(filename);
  b{i} = reshape(tmp, npf, length(tmp)/npf);
end

load('assemblelRHS.mat')
R = reshape(R, [npf length(R(:))/npf]);

for i = 1:4
  for j = 1:size(b{i},2)
    bj = b{i}(:,j);
    fj = dmd{i}.f2t(:,j);
    e1 = fj(1);
    l1 = fj(2);
    e2 = fj(3);
    l2 = fj(4);
    g1 = dmd{i}.elempart(e1);
    g2 = 0;
    if e2>0
      g2 = dmd{i}.elempart(e2);      
    end
    k = t2f(l1,g1);
    bk = R(:,k);
    if norm(bj-bk)>1e-8
      [i j e1 e2 g1 g2]
      error("RHS vectors do not match...");      
    end
  end
end

% ne1 = sum(elempartpts(1:2));
% nf = size(f2t,2);
% for i = 1:nf
%   if f2t(1,i)>ne1 
%     error("domain decomposition is incorrect because number of elements is out of bound");    
%   end
%   if f2t(3,i)>ne1 
%     
%   end
% end

% nfe = size(t2f,1);
% facerecv = 0*elemrecv;
% for i = 1:length(elemrecv)
%   e = elemrecv(i);  
%   for j = 1:nfe
%     if t2f(j,e) > 0
%       facerecv(i) = j;
%     end
%   end
% end
% 
% ne1 = sum(elempartpts(1:2));
% 
% facerecv = 0*elemrecv;
% facesend = 0*elemsend;
% nelemsend = length(elemsend);
% for i = 1:nelemsend
%   esend = elemsend(i);  
%   erecv = elemrecv(i);  
%   for j = 1:nfe
%     fsend = t2f(j,esend);
%     for k = 1:nfe
%       frecv = t2f(k,erecv);
%       if (fsend == frecv)
%         facesend(i) = j;
%         facerecv(i) = k;        
%       end
%     end
%   end
% end
% 
% 
% %     e1 = f2t(1,f);
% %     e2 = f2t(3,f);
% %     if e2>ne1
% %     end
% %     if t2f(j,e) > 0
% %       facesend(i) = j;
% %     end
