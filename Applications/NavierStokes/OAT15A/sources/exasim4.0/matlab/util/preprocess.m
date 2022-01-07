function [master,mesh,app] = preprocess(master,mesh,app)

hybrid = 'hdg';

% get dimensions
dim    = mesh.nd;
ngv    = master.ngv;
ngf    = master.ngf;
npv    = master.npv;
npf    = master.npf;
ne     = size(mesh.t,1);

%--------------- update MASTER structure ----------------%

% volume shape functions and their derivatives 
master.shapvgdotshapvl  = zeros(npv*npv,ngv,dim+1);      
for d=1:dim+1
    master.shapvt(:,:,d) = master.shapvl(:,:,d)';
    master.shapvg(:,:,d) = master.shapvl(:,:,d)*diag(master.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            master.shapvgdotshapvl((ii-1)*npv+jj,:,d) = master.shapvg(jj,:,d).*master.shapvl(ii,:,1);                    
        end
    end            
end

% face shape functions and their derivatives 
master.shapfgdotshapfc  = zeros(npf*npf,ngf,dim);      
for d=1:dim
    master.shapft(:,:,d) = master.shapfc(:,:,d)';
    master.shapfg(:,:,d) = master.shapfc(:,:,d)*diag(master.gwfc);
    for ii=1:npf
        for jj = 1:npf
            master.shapfgdotshapfc((ii-1)*npf+jj,:,d) = master.shapfg(jj,:,d).*master.shapfc(ii,:,1);                    
        end
    end            
end

mesh.permgeom = mesh.perm(:,:,1);
master.shapmv = master.shapvt;
master.shapmf = master.shapft;
master.shapmh = master.shapmf;
master.permgeom = mesh.perm(:,:,1);

% --------------- update MESH structure ----------------%
ngrsiz = 800;
ngr    = ceil(ne/ngrsiz);
ngrne  = round(ne/ngr);      
nk = 1:ngrne:ne;
nb = [nk(1:end); [nk(2:end)-1,ne]];
mesh.nb    = nb;

% Reorder faces
%[mesh.f,mesh.t2f,mesh.f2f,mesh.flev] = faceordering(mesh.f, mesh.t2f); 
[facecon,f,t2f] = mkconnectivity(mesh.p, mesh.t, master.porder, mesh.elemtype, mesh.bndexpr, []);
%figure(2);clf;plotface(mesh.p, mesh.t, f);
[mesh.f,mesh.t2f,mesh.facecon] = reorderface(f,t2f,facecon,app.boundarycondition,master.perm);
%figure(3);clf;plotface(mesh.p, mesh.t, mesh.f);

mesh.bf = reshape(mesh.f(abs(mesh.t2f'),end),[size(master.perm,2) size(mesh.t,1)]);
mesh.f2f = mkf2f(mesh.f, mesh.t2f);

if strcmp(hybrid,'edg')
    mesh.f(:,end+1) = 1;
    %fprintf('\n --- EDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
elseif strcmp(hybrid,'hdg')    
    mesh.f(:,end+1) = 0;
    %fprintf('\n --- HDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
elseif strcmp(hybrid,'hedg')
    a = mesh.f(:,end);
    i = (a>0);
    mesh.f(:,end+1) = 0;
    mesh.f(i,end)   = 1;
    %fprintf('\n --- H-EDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
elseif ~isempty(hybrid) 
    mesh.f = feval(hybrid,mesh);
else 
    error('Hybrid flag is not valid');
end

[mesh.elcon,mesh.nsiz] = elconnectivities(mesh);
mesh.f(:,end) = [];
