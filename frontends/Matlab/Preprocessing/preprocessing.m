function [app,mesh,master,dmd] = preprocessing(app,mesh)

coupledinterface = 0;
coupledcondition = 0;
coupledboundarycondition = 0;
if isfield(mesh, 'interfacecondition')
  coupledinterface = find(mesh.interfacecondition > 0);
  coupledcondition = mesh.interfacecondition(coupledinterface);  
  coupledboundarycondition = mesh.boundarycondition(coupledinterface);  
  if length(find(mesh.boundarycondition == coupledboundarycondition)) > 1
    error("mesh.boundarycondition is invalid because interface condition is not unique.");
  end
  if length(coupledinterface) > 1
    error("mesh.interfacecondition is invalid because there are more than one coupled interface.");
  end
  if isempty(coupledinterface)
    coupledinterface = 0;
    coupledcondition = 0;
    coupledboundarycondition = 0;
  end
end
app.coupledinterface = coupledinterface;
app.coupledcondition = coupledcondition;
app.coupledboundarycondition = coupledboundarycondition;

if app.modelnumber==0
    strn = "";
else
    strn = num2str(app.modelnumber);
end

if ~exist(char(app.exasimpath + "/build/model"), 'dir')
    mkdir(char(app.exasimpath + "/build/model"));
end
if  ~exist(char(app.buildpath + "/datain" + strn), 'dir')
    mkdir(char(app.buildpath + "/datain" + strn));
end
if  ~exist(char(app.buildpath + "/dataout" + strn), 'dir')
    mkdir(char(app.buildpath + "/dataout" + strn));
end
filename = app.buildpath + "/datain" + strn + "/";
fileapp = filename + "app.bin";
filemaster = filename + "master.bin";
endian = 'native';

if app.preprocessmode==0    
    % update app structure    
    writeapp(app,fileapp,endian);       
    master = [];
    dmd = [];
    return;
end

app.nd  = size(mesh.p,1);
app.ncx = app.nd;
[app.nve,app.ne] = size(mesh.t);
app.elemtype = 0;
if (app.nd==2) && (app.nve==4)
    app.elemtype=1;    
end
if (app.nd==3) && (app.nve==8)
    app.elemtype=1;    
end

if isfield(app, 'pgauss')==0
  app.pgauss = 2*app.porder;
elseif app.pgauss < 2*app.porder
  app.pgauss = 2*app.porder;
end

% master struct
master = Master(app);
writemaster(master,filemaster,'native');        

mesh.coupledinterface = app.coupledinterface;
mesh.porder = app.porder;
mesh.perm = master.perm;

% obtain the PDE model
pdemodel = str2func(app.modelfile);
pde = pdemodel();

app.boundaryconditions = mesh.boundarycondition;
app.uinf = app.externalparam;
nuinf = length(app.uinf);
nparam = length(app.physicsparam);
xdgsym = sym('xdg',[app.ncx 1]); 
uinfsym = sym('uinf',[nuinf 1]); 
paramsym = sym('param',[nparam 1]); 
if isfield(pde, 'initu')    
    udgsym = pde.initu(xdgsym, paramsym, uinfsym);         
    app.ncu = length(udgsym(:));
else
    error("pde.initu is not defined");
end
if isfield(pde, 'initv')
    odgsym = pde.initv(xdgsym, paramsym, uinfsym);         
    app.nco = length(odgsym(:));
elseif isfield(mesh, 'vdg')            
    app.nco = size(mesh.vdg,2);    
else    
    app.nco = 0;
end
if isfield(pde, 'initw')
    wdgsym = pde.initw(xdgsym, paramsym, uinfsym);         
    app.ncw = length(wdgsym(:));
elseif isfield(mesh, 'wdg')            
    app.ncw = size(mesh.wdg,2);        
else    
    app.ncw = 0;
end

if app.model=="ModelC" || app.model=="modelC"
    app.wave = 0;
    app.nc = app.ncu;
elseif app.model=="ModelD" || app.model=="modelD"     
    app.wave = 0;
    app.nc = (app.ncu)*(app.nd+1);
elseif app.model=="ModelW" || app.model=="modelW"
    app.tdep = 1;
    app.wave = 1;
    app.nc = (app.ncu)*(app.nd+1);
end
app.ncq = app.nc - app.ncu;
app.nch  = app.ncu;                

if max(app.dt)>0
    app.tdep = 1;
else
    app.tdep = 0;
end

disp('run facenumbering...');  
[mesh.f, mesh.tprd, t2t] = facenumbering(mesh.p,mesh.t,app.elemtype,mesh.boundaryexpr,mesh.periodicexpr);

mpiprocs = app.mpiprocs;
%dmd = meshpartition(mesh.p,mesh.t,mesh.f,t2t,mesh.tprd,app.elemtype,app.boundaryconditions,mesh.boundaryexpr,mesh.periodicexpr,app.porder,mpiprocs,app.metis);
if (app.hybrid ==1)
  dmd = meshpartitionhdg(mesh.tprd,mesh.f,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,coupledinterface,mpiprocs,app.metis);
else
  dmd = meshpartition2(mesh.tprd,mesh.f,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,mpiprocs,app.metis);
end

for i = 1:mpiprocs       
    % create DG nodes
    if isfield(mesh, 'dgnodes')
        xdg = mesh.dgnodes(:,:,dmd{i}.elempart);
    else
        xdg = createdgnodes(mesh.p,mesh.t(:,dmd{i}.elempart),mesh.f(:,dmd{i}.elempart),mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);    
        if mpiprocs==1
          mesh.dgnodes = xdg;        
        end
    end    
        
    %xdg = createdgnodes(mesh.p,mesh.t(:,dmd{i}.elempart),mesh.f(:,dmd{i}.elempart),mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);    
    
    [~,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(xdg,1e-6);
    disp(['Writing initial solution into file ' num2str(i) '...']);
    if (mpiprocs>1) || (coupledinterface>0)
        fileID1 = fopen(filename + "sol" + string(i) + ".bin",'w');
    else
        fileID1 = fopen(filename + "sol" + ".bin",'w');
    end
    ndims = zeros(12,1);           
    ndims(1) = length(dmd{i}.elempart); % number of elements
    ndims(2) = sum(dmd{i}.facepartpts); % number of faces
    ndims(3) = size(master.perm,2); % number of faces per element          
    ndims(4) = master.npe;
    ndims(5) = master.npf;            
    ndims(6) = app.nc;
    ndims(7) = app.ncu;
    ndims(8) = app.ncq;
    ndims(9) = app.ncw;
    ndims(10) = app.nco;
    ndims(11) = app.nch;
    ndims(12) = app.ncx;

    nsize = zeros(20,1);
    nsize(1) = length(ndims(:));
    nsize(2) = length(xdg(:));
    %nsize(3) = length(udg(:)); 
%     nsize(4) = length(odg(:));    
%     nsize(5) = length(wdg(:));    
    if isfield(mesh, 'udg')        
        nsize(3) = numel(mesh.udg(:,:,dmd{i}.elempart));
    end
    if isfield(mesh, 'vdg')        
        nsize(4) = numel(mesh.vdg(:,:,dmd{i}.elempart));
    end
    if isfield(mesh, 'wdg')        
        nsize(5) = numel(mesh.wdg(:,:,dmd{i}.elempart));
    end
    if isfield(app, 'read_uh') 
        if app.read_uh
            % First, read from UH (specifically np{i})
            % Just write that directly, no need for UH(i,j,k);
            if isfield(app, 'read_uh_steady')
                disp("FORCE THE READ OF UH STEADY")
                fileID = fopen(app.buildpath+"/dataout/out_uhat_np"+string(i-1)+".bin",'r');
            else
            fileID = fopen(app.buildpath+"/dataout/out_uhat_t1500_np"+string(i-1)+".bin",'r');
            end
            UH_tmp = fread(fileID,'double');
            nsize(6) = length(UH_tmp); 
        end
    end

    fwrite(fileID1,length(nsize(:)),'double',endian);
    fwrite(fileID1,nsize(:),'double',endian);
    fwrite(fileID1,ndims(:),'double',endian);
    fwrite(fileID1,xdg(:),'double',endian);
    %fwrite(fileID1,udg(:),'double',endian);        
%     fwrite(fileID1,odg(:),'double',endian);
%     fwrite(fileID1,wdg(:),'double',endian);
    if isfield(mesh, 'udg')        
        fwrite(fileID1,mesh.udg(:,:,dmd{i}.elempart),'double',endian);                
    end
    if isfield(mesh, 'vdg')        
        fwrite(fileID1,mesh.vdg(:,:,dmd{i}.elempart),'double',endian);                
    end
    if isfield(mesh, 'wdg')        
        fwrite(fileID1,mesh.wdg(:,:,dmd{i}.elempart),'double',endian);                
    end
    if isfield(app, 'read_uh')    
        if app.read_uh 
        % Just write that directly, no need for UH(i,j,k);
        fwrite(fileID1,UH_tmp,'double',endian);        
        end
    end
    fclose(fileID1);         
    
    % divide elements and faces into blocks
    if mpiprocs==1
        if coupledinterface>0
          me = cumsum([0 dmd{i}.intepartpts(1) dmd{i}.intepartpts(2)]);
          [eblks,nbe] = mkfaceblocks(me,[0 -1],app.neb);                
        else
          ne = length(dmd{i}.elempart);        
          [eblks,nbe] = mkelemblocks(ne,app.neb);
          eblks(3,:) = 0;        
        end     
        mf = cumsum([0 dmd{i}.facepartpts]);    
        [fblks,nbf] = mkfaceblocks(mf,dmd{i}.facepartbnd,app.nfb);       
        neb = max(eblks(2,:)-eblks(1,:))+1;
        nfb = max(fblks(2,:)-fblks(1,:))+1;        
    else
        if coupledinterface>0
          me = cumsum([0 dmd{i}.intepartpts(1) dmd{i}.intepartpts(2) dmd{i}.intepartpts(3) dmd{i}.intepartpts(4)]);
          [eblks,nbe] = mkfaceblocks(me,[0 -1 1 2],app.neb);          
        else
          me = cumsum([0 dmd{i}.elempartpts(1) dmd{i}.elempartpts(2) dmd{i}.elempartpts(3)]);
          [eblks,nbe] = mkfaceblocks(me,[0 1 2],app.neb);          
        end
        mf = cumsum([0 dmd{i}.facepartpts]);                 
        [fblks,nbf] = mkfaceblocks(mf,dmd{i}.facepartbnd,app.nfb);        
        neb = max(eblks(2,:)-eblks(1,:))+1;
        nfb = max(fblks(2,:)-fblks(1,:))+1;        
    end            
    npe = master.npe;
    nfe = size(master.perm,2);
    facecon1 = reshape(dmd{i}.facecon(:,1,:),[size(dmd{i}.facecon,1) size(dmd{i}.facecon,3)]);
    facecon2 = reshape(dmd{i}.facecon(:,2,:),[size(dmd{i}.facecon,1) size(dmd{i}.facecon,3)]);      
    ind = [];        
    for ii = 1:size(fblks,2)
        if fblks(3,ii)>0
            ind = [ind fblks(1,ii):fblks(2,ii)];
        end
    end            
    facecon2(:,ind)=[];        
    [rowe2f1,cole2f1,ent2ind1] = mkdge2dgf(facecon1,npe*length(dmd{i}.elempart));                
    [rowe2f2,cole2f2,ent2ind2] = mkdge2dgf(facecon2,npe*length(dmd{i}.elempart));      
    
    %writebin("facecon" + num2str(i) + ".bin",facecon2(:));
        
    disp(['Writing mesh into file ' num2str(i) '...']); 
    if (mpiprocs>1) || (coupledinterface>0)
        fileID2 = fopen(filename + "mesh" + string(i) + ".bin",'w');
    else
        fileID2 = fopen(filename + "mesh" + ".bin",'w');
    end
    ndims = zeros(20,1);
    ndims(1) = size(mesh.p,1);
    ndims(2) = length(dmd{i}.elempart);
    ndims(3) = sum(dmd{i}.facepartpts);
    ndims(4) = max(mesh.t(:));
    ndims(5) = nfe;
    ndims(6) = nbe;
    ndims(7) = neb;
    ndims(8) = nbf;
    ndims(9) = nfb;
        
    nsize = zeros(50,1);
    nsize(1) = length(ndims(:));
    nsize(2) = length(dmd{i}.facecon(:));  
    nsize(3) = length(eblks(:)); 
    nsize(4) = length(fblks(:)); 
    nsize(5) = length(dmd{i}.nbsd(:)); 
    nsize(6) = length(dmd{i}.elemsend(:)); 
    nsize(7) = length(dmd{i}.elemrecv(:)); 
    nsize(8) = length(dmd{i}.elemsendpts(:)); 
    nsize(9) = length(dmd{i}.elemrecvpts(:));             
    nsize(10) = length(dmd{i}.elempart(:)); 
    nsize(11) = length(dmd{i}.elempartpts(:)); 
    nsize(12) = length(cgelcon(:));  
    nsize(13) = length(rowent2elem(:));  
    nsize(14) = length(cgent2dgent(:));  
    nsize(15) = length(colent2elem(:));                          
    nsize(16) = length(rowe2f1(:));  
    nsize(17) = length(cole2f1(:));  
    nsize(18) = length(ent2ind1(:));                          
    nsize(19) = length(rowe2f2(:));  
    nsize(20) = length(cole2f2(:));  
    nsize(21) = length(ent2ind2(:));  
    if (app.hybrid > 0)
      nsize(22) = length(dmd{i}.f2t(:));  
      nsize(23) = length(dmd{i}.elemcon(:));  
      nsize(24) = length(master.perm(:));  
      nsize(25) = length(dmd{i}.bf(:));           
%       if coupledinterface>0
%         nsize(26) = length(dmd{i}.intepartpts(:));    
%       end
    end
    fwrite(fileID2,length(nsize(:)),'double',endian);
    fwrite(fileID2,nsize(:),'double',endian);
    fwrite(fileID2,ndims(:),'double',endian);
    fwrite(fileID2,reshape(permute(dmd{i}.facecon,[2 1 3]),[2*master.npf(1)*ndims(3) 1]),'double',endian);
    fwrite(fileID2,eblks(:),'double',endian);
    fwrite(fileID2,fblks(:),'double',endian);
    fwrite(fileID2,dmd{i}.nbsd(:),'double',endian);
    fwrite(fileID2,dmd{i}.elemsend(:),'double',endian);
    fwrite(fileID2,dmd{i}.elemrecv(:),'double',endian);
    fwrite(fileID2,dmd{i}.elemsendpts(:),'double',endian);
    fwrite(fileID2,dmd{i}.elemrecvpts(:),'double',endian);
    fwrite(fileID2,dmd{i}.elempart(:),'double',endian);
    fwrite(fileID2,dmd{i}.elempartpts(:),'double',endian);
    fwrite(fileID2,cgelcon(:)-1,'double',endian);
    fwrite(fileID2,rowent2elem(:),'double',endian);
    fwrite(fileID2,cgent2dgent(:)-1,'double',endian);
    fwrite(fileID2,colent2elem(:)-1,'double',endian);            
    fwrite(fileID2,rowe2f1(:),'double',endian);
    fwrite(fileID2,cole2f1(:)-1,'double',endian);
    fwrite(fileID2,ent2ind1(:)-1,'double',endian);                 
    fwrite(fileID2,rowe2f2(:),'double',endian);
    fwrite(fileID2,cole2f2(:)-1,'double',endian);
    fwrite(fileID2,ent2ind2(:)-1,'double',endian);        
    if (app.hybrid > 0) % HDG
      fwrite(fileID2,dmd{i}.f2t(:)-1,'double',endian);        
      fwrite(fileID2,dmd{i}.elemcon(:)-1,'double',endian);        
      fwrite(fileID2,master.perm(:)-1,'double',endian);        
      fwrite(fileID2,dmd{i}.bf(:),'double',endian);      
%       if coupledinterface>0
%         fwrite(fileID2,dmd{i}.intepartpts(:),'double',endian);      
%       end
    end
    fclose(fileID2);             
end

app = writeapp(app,fileapp,endian);                

mesh.telem = master.telem;
mesh.tface = master.telem;
mesh.xpe = master.xpe;
mesh.xpf = master.xpf;
% for i = 1:mpiprocs   
%     dmd{i}.elem2cpu = [];
%     dmd{i}.elemrecv = [];
%     dmd{i}.nbsd = [];
%     dmd{i}.elemsend = [];
%     dmd{i}.elemsendpts = [];
%     dmd{i}.elemrecvpts = [];
%     dmd{i}.facepartpts = [];
%     dmd{i}.facepartbnd = [];
%     dmd{i}.facecon = [];
% end
