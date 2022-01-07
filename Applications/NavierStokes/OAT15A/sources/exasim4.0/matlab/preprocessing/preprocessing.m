function [app,dmd] = preprocessing(app,p,t,dgnodes,UDG,ODG,WDG,bndexpr,periodicexpr,endianType,preprocessmode)

if nargin < 11; preprocessmode = 1; end
if nargin < 10; endianType = 0; end
if endianType == 0; endian = 'native';
elseif endianType == 1; endian = 'ieee-le';
elseif endianType == 2; endian = 'ieee-be';
end

app.nd  = size(p,2);
app.nc  = size(UDG,2);          
app.ncx = size(dgnodes,2);
app.nco = size(ODG,2);
app.ncp = 0;

if app.pdemodel==1
    app.wave = 0;
    app.ncu = app.nc;
elseif app.pdemodel==2    
    app.wave = 0;
    app.ncu = (app.nc-app.ncp)/(app.nd+1);
elseif app.pdemodel==3
    app.tdep = 1;
    app.wave = 1;
    app.ncu = (app.nc-app.ncp)/(app.nd+1);
end

app.ncq = app.nc - app.ncu - app.ncp;
app.nch  = app.ncu;                

if max(app.dt)>0
    app.tdep = 1;
else
    app.tdep = 0;
end

app.pgauss = 2*app.porder;
app.neb = 512*32;  % number of elements per block
app.nfb = 512*64; % number of faces per block
%app

if isfield(app,'filedir') == 0
    filename = app.appname;
else
    filename = app.filedir;
    mkdir(filename);
end

mpiprocs = app.mpiprocs;
if preprocessmode==0    
    % update app structure
    fileapp = [filename 'app.bin'];
    app = writeapp(app,fileapp,endian);        
else        
    % master strcuture
    filemaster = [filename 'master.bin'];
    fileID = fopen(filemaster,'w');
    master = mkmasterelement(app.nd,app.porder,app.porder,app.pgauss,app.pgauss,app.elemtype,app.nodetype);        
    writemaster(master,fileID,endian);        
    fclose(fileID);
    
    disp('run mkconnectivity...');  
    [facecon,f,t2f,t2t,t] = mkconnectivity(p, t, master.porder, master.elemtype, bndexpr, periodicexpr);
    %save temp.mat f t2f facecon t t2t;
        
    npe = size(UDG,1);
    disp('run reorderface...');    
    if mpiprocs==1       
    [f,t2f,facecon,mf,boundarycondition] = reorderface(f,t2f,facecon,app.boundarycondition,master.perm);          
    else
    [f,t2f,facecon,mf,boundarycondition] = reorderface_mpi(f,t2f,facecon,app.boundarycondition,master.perm);            
    end
    %save temp.mat f t2f facecon mf;
    app.ne = size(t2f,1);
    app.nf = size(f,1);
    app.nv = size(p,1);        
    if mpiprocs>1        
        disp('run mkpartition...');  
        dmd = mkpartition(facecon,f,t2f,t2t,t,mpiprocs,npe,app.boundarycondition,master.perm);        
%         save test.mat dmd;        
%         plotpartition(p,t,f,dmd);
        
        for i = 1:mpiprocs            
            xdg = dgnodes(:,:,dmd{i}.elempart);                        
            if isempty(UDG)==0
                udg = UDG(:,:,dmd{i}.elempart);
            else
                udg = [];
            end
            if isempty(ODG)==0
                odg = ODG(:,:,dmd{i}.elempart);
            else
                odg = [];
            end
            if isempty(WDG)==0
                wdg = WDG(:,:,dmd{i}.elempart);
            else
                wdg = [];
            end
            
            [~,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(xdg(:,1:app.nd,:),1e-6);

            disp(['Writing initial solution into file ' num2str(i) '...']); disp(' ');
            fileID1 = fopen([filename,'sol' num2str(i) '.bin'],'w');
            ndims = zeros(12,1);           
            ndims(1) = length(dmd{i}.elempart); % number of elements
            ndims(2) = length(dmd{i}.facepart); % number of faces
            ndims(3) = size(t2f,2); % number of faces per element          
            ndims(4) = master.npv;
            ndims(5) = master.npf(1);            
            ndims(6) = app.nc;
            ndims(7) = app.ncu;
            ndims(8) = app.ncq;
            ndims(9) = app.ncp;
            ndims(10) = app.nco;
            ndims(11) = app.nch;
            ndims(12) = app.ncx;
            
            nsize = zeros(20,1);
            nsize(1) = length(ndims(:));
            nsize(2) = length(xdg(:));
            nsize(3) = length(udg(:)); 
            nsize(4) = length(odg(:));    
            nsize(5) = length(wdg(:));    
            
            fwrite(fileID1,length(nsize(:)),'double',endian);
            fwrite(fileID1,nsize(:),'double',endian);
            fwrite(fileID1,ndims(:),'double',endian);
            fwrite(fileID1,xdg(:),'double',endian);
            fwrite(fileID1,udg(:),'double',endian);        
            fwrite(fileID1,odg(:),'double',endian);
            fwrite(fileID1,wdg(:),'double',endian);
            fclose(fileID1);         
            
            % divide elements and faces into blocks
            me = cumsum([0 dmd{i}.elempartpts(1) dmd{i}.elempartpts(2) dmd{i}.elempartpts(3)]);
            [eblks,nbe] = mkfaceblocks(me,[0 1 2],app.neb);                           
            mf = cumsum([0 dmd{i}.facepartpts]);                 
            [fblks,nbf] = mkfaceblocks(mf,dmd{i}.facepartbnd,app.nfb);        
            neb = max(eblks(2,:)-eblks(1,:))+1;
            nfb = max(fblks(2,:)-fblks(1,:))+1;            
            %nbf0 = find(fblks(2,:)==dmd{i}.nfbind(1));
            %nbf1 = find(fblks(2,:)==dmd{i}.nfbind(2));   
            eblksa{i} = eblks;
            fblksa{i} = fblks;
            save dmdblks.mat eblksa fblksa
            
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
            
            disp(['Writing mesh into file ' num2str(i) '...']); disp(' ');
            fileID2 = fopen([filename,'mesh' num2str(i) '.bin'],'w');    
            ndims = zeros(20,1);
            ndims(1) = size(p,2);
            ndims(2) = length(dmd{i}.elempart);
            ndims(3) = length(dmd{i}.facepart);
            ndims(4) = max(max(t(dmd{i}.elempart,:)));
            ndims(5) = size(t2f,2);
            ndims(6) = nbe;
            ndims(7) = neb;
            ndims(8) = nbf;
            ndims(9) = nfb;
            %ndims(10) = nbf0;
            %ndims(11) = nbf1;

            nsize = zeros(30,1);
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
            %nsize(12) = length(dmd{i}.facecon(:));  
            %nsize(13) = length(dmd{i}.t2f(:));  
            %nsize(14) = length(dmd{i}.f(:));  
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
            %[~,cgelcon,rowent2elem,~,cgent2dgent]
            %[rowe2f1,cole2f1,ent2ind1] = mkdge2dgf(facecon1);
            %[rowe2f2,cole2f2,ent2ind2] = mkdge2dgf(facecon2);
            
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
            %fwrite(fileID2,dmd{i}.facecon(:),'double',endian);
            %fwrite(fileID2,dmd{i}.t2f(:),'double',endian);
            %fwrite(fileID2,dmd{i}.f(:),'double',endian);
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
            fclose(fileID2);         
        end
    else        
        dmd = 0;
        disp('Writing initial solution into file ...'); disp(' ');
        fileID = fopen([filename,'sol.bin'],'w');
        
        ndims = zeros(12,1);           
        ndims(1) = size(t,1); % number of elements
        ndims(2) = size(f,1); % number of faces
        ndims(3) = size(t2f,2); % number of faces per element          
        ndims(4) = master.npv;
        ndims(5) = master.npf(1);            
        ndims(6) = app.nc;
        ndims(7) = app.ncu;
        ndims(8) = app.ncq;
        ndims(9) = app.ncp;
        ndims(10) = app.nco;
        ndims(11) = app.nch;
        ndims(12) = app.ncx;
        
        nsize = zeros(20,1);
        nsize(1) = length(ndims(:));
        nsize(2) = length(dgnodes(:));
        nsize(3) = length(UDG(:)); 
        nsize(4) = length(ODG(:));         
        nsize(5) = length(WDG(:));         
            
        fwrite(fileID,length(nsize(:)),'double',endian);
        fwrite(fileID,nsize(:),'double',endian);
        fwrite(fileID,ndims(:),'double',endian);    
        fwrite(fileID,dgnodes(:),'double',endian);
        fwrite(fileID,UDG(:),'double',endian);        
        fwrite(fileID,ODG(:),'double',endian);   
        fwrite(fileID,WDG(:),'double',endian);   
        fclose(fileID);      
        
        [~,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(dgnodes(:,1:app.nd,:),1e-5);
        
        % divide elements and faces into blocks
        ne = size(t2f,1);
        [eblks,nbe] = mkelemblocks(ne,app.neb);
        eblks(3,:) = 0;        
        [fblks,nbf] = mkfaceblocks(mf,[0 boundarycondition(:)'],app.nfb);       
        neb = max(eblks(2,:)-eblks(1,:))+1;
        nfb = max(fblks(2,:)-fblks(1,:))+1;
                        
        facecon1 = reshape(facecon(:,1,:),[size(facecon,1) size(facecon,3)]);
        facecon2 = reshape(facecon(:,2,:),[size(facecon,1) size(facecon,3)]);      
        ind = [];        
        for i = 1:size(fblks,2)
            if fblks(3,i)>0
                ind = [ind fblks(1,i):fblks(2,i)];
            end
        end        
        facecon2(:,ind)=[];        
        [rowe2f1,cole2f1,ent2ind1] = mkdge2dgf(facecon1,npe*ne);                
        [rowe2f2,cole2f2,ent2ind2] = mkdge2dgf(facecon2,npe*ne);        
        
        disp('Writing mesh into file...'); disp(' ');
        fileID = fopen([filename,'mesh.bin'],'w');    
        ndims = zeros(20,1);
        ndims(1) = size(p,2);
        ndims(2) = size(t,1);
        ndims(3) = size(f,1);
        ndims(4) = max(t(:));
        ndims(5) = size(t2f,2);
        ndims(6) = nbe;
        ndims(7) = neb;
        ndims(8) = nbf;
        ndims(9) = nfb;
        
        nsize = zeros(20,1);
        nsize(1) = length(ndims(:));
        nsize(2) = length(facecon(:));  
        nsize(3) = length(eblks(:)); 
        nsize(4) = length(fblks(:)); 
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

        % write mesh structure to files
        fwrite(fileID,length(nsize(:)),'double',endian);
        fwrite(fileID,nsize(:),'double',endian);
        fwrite(fileID,ndims(:),'double',endian);
        fwrite(fileID,reshape(permute(facecon,[2 1 3]),[2*master.npf(1)*size(f,1) 1]),'double',endian);
        fwrite(fileID,eblks(:),'double',endian);
        fwrite(fileID,fblks(:),'double',endian);        
        fwrite(fileID,cgelcon(:)-1,'double',endian);
        fwrite(fileID,rowent2elem(:),'double',endian);
        fwrite(fileID,cgent2dgent(:)-1,'double',endian);
        fwrite(fileID,colent2elem(:)-1,'double',endian);   
        fwrite(fileID,rowe2f1(:),'double',endian);
        fwrite(fileID,cole2f1(:)-1,'double',endian);
        fwrite(fileID,ent2ind1(:)-1,'double',endian);                 
        fwrite(fileID,rowe2f2(:),'double',endian);
        fwrite(fileID,cole2f2(:)-1,'double',endian);
        fwrite(fileID,ent2ind2(:)-1,'double',endian);                 
        fclose(fileID);                    
    end            

    fileapp = [filename 'app.bin'];
    app = writeapp(app,fileapp,endian);                
end


function app = writeapp(app,filename,endian)

if strcmp(app.appname,'euler')
    appname = 0;    
elseif strcmp(app.appname,'ns')
    appname = 1;        
elseif strcmp(app.appname,'poisson')
    appname = 2;   
elseif strcmp(app.appname,'ransSA')
    appname = 3;  
elseif strcmp(app.appname,'leuq')
    appname = 4;      
elseif strcmp(app.appname,'ledisp')
    appname = 5;
elseif strcmp(app.appname,'leuqNH')
    appname = 6;
elseif strcmp(app.appname,'generic')
    appname = 1000;    
else
    appname = 1000;
end
if isfield(app,'hybrid') == 0
    hybridn = 0;
else
    hybridn = app.hybrid;
end
if isfield(app,'wave') == 0
    app.wave = 0;
end
if isfield(app,'linearproblem') == 0
    app.linearproblem = 0;
end
if isfield(app,'debugmode') == 0
    app.debugmode = 0;
end
if isfield(app,'temporalscheme') == 0
    app.temporalscheme = 0;
end
if isfield(app,'torder') == 0
    app.torder = 1;
end
if isfield(app,'nstage') == 0
    app.nstage = 1;
end
if isfield(app,'convStabMethod') == 0
    app.convStabMethod = 0;
end
if isfield(app,'diffStabMethod') == 0
    app.diffStabMethod = 0;
end
if isfield(app,'rotatingFrame') == 0
    app.rotatingFrame = 0;
end
if isfield(app,'viscosityModel') == 0
    app.viscosityModel = 0;
end
if isfield(app,'SGSmodel') == 0
    app.SGSmodel = 0;
end
if isfield(app,'ALE') == 0
    app.ALE = 0;
end
if isfield(app,'AV') == 0
    app.AV = 0;
end
if isfield(app,'nonlinearsolver') == 0
    app.nonlinearsolver = 0;
end
if isfield(app,'linearsolver') == 0
    app.linearsolver = 0;
end
if isfield(app,'PTCiter') == 0
    app.PTCiter = 20;
end
if isfield(app,'PTCtol') == 0
    app.PTCtol = 1e-6;
end
if isfield(app,'linearsolveriter') == 0
    app.linearsolveriter = 200;
end
if isfield(app,'GMRESrestart') == 0
    app.GMRESrestart = 25;
end
if isfield(app,'linearsolvertol') == 0
    app.linearsolvertol = 1e-3;
end
if isfield(app,'GMRESortho') == 0
    app.GMRESortho = 0;
end
if isfield(app,'preconditioner') == 0
    app.preconditioner = 0;
end
if isfield(app,'precMatrixType') == 0
    app.precMatrixType = 0;
end
if isfield(app,'ptcMatrixType') == 0
    app.ptcMatrixType = 0;
end
if isfield(app,'runmode') == 0
    app.runmode = 0;
end
if isfield(app,'tdfunc') == 0
    app.tdfunc = 1;
end
if isfield(app,'source') == 0
    app.source = 1;
end

if isfield(app,'matvectol') == 0
    app.matvectol = 1e-6;
end
if isfield(app,'matvecorder') == 0
    app.matvecorder = 1;
end
if isfield(app,'PTCparam') == 0
    app.PTCparam = 1;
end
if isfield(app,'RBdim') == 0
    app.RBdim = 5;
end
if isfield(app,'saveSolFreq') == 0
    app.saveSolFreq = 1;
end
if isfield(app,'saveSolOpt') == 0
    app.saveSolOpt = 1;
end

% if isfield(app,'dtcoef_u') == 0
%     app.dtcoef_u = 0;
% end
% if isfield(app,'dtcoef_q') == 0
%     app.dtcoef_q = 0;
% end
% if isfield(app,'dtcoef_p') == 0
%     app.dtcoef_p = 0;
% end
if isfield(app,'time') == 0  || isempty(app.time)==1
    app.time = 0;
end
if isfield(app,'nfile') == 0
    app.nfile = 1;      
end
if isfield(app,'flag') == 0
    app.flag = [];         
end
if isfield(app,'factor') == 0
    app.factor = [];         
end
if isfield(app,'problem') == 0
    app.problem = [];         
end
if isfield(app,'solversparam') == 0
    app.solversparam = [];         
end
if isfield(app,'tau') == 0
    app.tau = [];         
end
if isfield(app,'tdepbc') == 0
    app.tdepbc = [];         
else
    app.tdepbc = unique(app.tdepbc(:));
end
% if length(app.dtcoef_u)==1
%     app.dtcoef_u = app.dtcoef_u*ones(1,app.ncu);
% end
% if app.ncq>0
%     if length(app.dtcoef_q)==1
%         app.dtcoef_q = app.dtcoef_q*ones(1,app.ncq);
%     end
%     if max(app.dtcoef_q)<=0
%         error('dtcoef_q must be positive');
%     end
% else
%     app.dtcoef_q = [];
% end
% if app.ncp>0
%     if length(app.dtcoef_p)==1
%         app.dtcoef_p = app.dtcoef_p*ones(1,app.ncp);
%     end
%     if max(app.dtcoef_p)<=0
%         error('dtcoef_p must be positive');
%     end
% else
%     app.dtcoef_p = [];
% end
% if app.tdep==1
%     if max(app.dtcoef_u)<=0
%         error('dtcoef_u must be positive for time-dependent problems');
%     end
% end

app.flag   = [app.tdep app.wave app.linearproblem app.debugmode app.matvecorder app.GMRESortho...  
              app.preconditioner app.precMatrixType app.ptcMatrixType app.runmode app.tdfunc app.source app.flag];
app.problem  = [hybridn appname app.temporalscheme app.torder app.nstage app.convStabMethod...
               app.diffStabMethod app.rotatingFrame app.viscosityModel app.SGSmodel app.ALE app.AV...
               app.linearsolver app.PTCiter app.linearsolveriter app.GMRESrestart app.RBdim ...
               app.saveSolFreq app.saveSolOpt app.problem];
%app.factor = [app.dtcoef_u app.dtcoef_q app.dtcoef_p app.time app.factor];           
app.factor = [app.time app.factor];           
%app.physicsparam  = reshape(cell2mat(app.arg),[],1); 
app.solversparam = [app.PTCtol app.linearsolvertol app.matvectol app.PTCparam app.solversparam];        

ndims = zeros(40,1);
ndims(1) = app.mpiprocs;  % number of processors
ndims(2) = app.nd;
ndims(3) = 0;
ndims(4) = 0;
ndims(5) = 0;
ndims(6) = app.nc;
ndims(7) = app.ncu;
ndims(8) = app.ncq;
ndims(9) = app.ncp;
ndims(10) = app.nco;
ndims(11) = app.nch;
ndims(12) = app.ncx;

nsize = zeros(10,1);
nsize(1) = length(ndims(:));
nsize(2) = length(app.flag(:));  % length of flag
nsize(3) = length(app.problem(:)); % length of physics
nsize(4) = length(app.uinf(:)); % boundary data
nsize(5) = length(app.dt(:)); % number of time steps
nsize(6) = length(app.factor(:)); % length of factor
nsize(7) = length(app.physicsparam(:)); % number of physical parameters
nsize(8) = length(app.solversparam(:)); % number of solver parameters
nsize(9) = length(app.tau(:)); % number of stabilization parameters
nsize(10) = length(app.tdepbc(:)); % number of time-dependent boundary conditions

app.nsize = nsize;
app.ndims = ndims;
disp('Writing app into file...'); disp(' ');
fileID = fopen(filename,'w');
fwrite(fileID,length(app.nsize(:)),'double',endian);
fwrite(fileID,app.nsize(:),'double',endian);
fwrite(fileID,app.ndims(:),'double',endian);
fwrite(fileID,app.flag(:),'double',endian);
fwrite(fileID,app.problem(:),'double',endian);
fwrite(fileID,app.uinf(:),'double',endian);
fwrite(fileID,app.dt(:),'double',endian);
fwrite(fileID,app.factor(:),'double',endian);
fwrite(fileID,app.physicsparam(:),'double',endian);
fwrite(fileID,app.solversparam(:),'double',endian);
fwrite(fileID,app.tau(:),'double',endian);
fwrite(fileID,app.tdepbc(:),'double',endian);
fclose(fileID);

function writemaster(master,fileID,endian)

disp('Writing master into file...'); disp(' ');

ndims = zeros(20,1);
ndims(1) = master.nd;
ndims(2) = master.elemtype;
ndims(3) = master.nodetype;
ndims(4) = master.porder;
ndims(5) = master.pgauss;
ndims(6) = master.npv;
ndims(7) = master.npf(1);
ndims(8) = master.ngv;
ndims(9) = master.ngf(1);
ndims(10) = length(master.ploc1d);
ndims(11) = length(master.gp1d);

nsize = zeros(22,1);
nsize(1) = length(ndims(:));
nsize(2) = length(master.shapvt(:));  
nsize(3) = length(master.shapvg(:)); 
nsize(4) = length(master.shapft{1}(:)); 
nsize(5) = length(master.shapfg{1}(:));
nsize(6) = length(master.shapnt(:));  
nsize(7) = length(master.shapnl(:)); 
nsize(8) = length(master.shapfnt{1}(:)); 
nsize(9) = length(master.shapfn{1}(:));
nsize(10) = length(master.plocvl(:)); 
nsize(11) = length(master.gpvl(:)); 
nsize(12) = length(master.gwvl(:)); 
nsize(13) = length(master.plocfc{1}(:)); 
nsize(14) = length(master.gpfc{1}(:)); 
nsize(15) = length(master.gwfc{1}(:)); 
nsize(16) = length(master.shap1dgt(:)); 
nsize(17) = length(master.shap1dgw(:)); 
nsize(18) = length(master.shap1dnt(:)); 
nsize(19) = length(master.shap1dnl(:));
nsize(20) = length(master.ploc1d(:)); 
nsize(21) = length(master.gp1d(:)); 
nsize(22) = length(master.gw1d(:)); 


% write master structure to files
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,ndims(:),'double',endian);
fwrite(fileID,master.shapvt(:),'double',endian);
fwrite(fileID,master.shapvg(:),'double',endian);
fwrite(fileID,master.shapft{1}(:),'double',endian);
fwrite(fileID,master.shapfg{1}(:),'double',endian);
fwrite(fileID,master.shapnt(:),'double',endian);
fwrite(fileID,master.shapnl(:),'double',endian);
fwrite(fileID,master.shapfnt{1}(:),'double',endian);
fwrite(fileID,master.shapfn{1}(:),'double',endian);
fwrite(fileID,master.plocvl(:),'double',endian);
fwrite(fileID,master.gpvl(:),'double',endian);
fwrite(fileID,master.gwvl(:),'double',endian);
fwrite(fileID,master.plocfc{1}(:),'double',endian);
fwrite(fileID,master.gpfc{1}(:),'double',endian);
fwrite(fileID,master.gwfc{1}(:),'double',endian);
fwrite(fileID,master.shap1dgt(:),'double',endian);
fwrite(fileID,master.shap1dgw(:),'double',endian);
fwrite(fileID,master.shap1dnt(:),'double',endian);
fwrite(fileID,master.shap1dnl(:),'double',endian);
fwrite(fileID,master.ploc1d(:),'double',endian);
fwrite(fileID,master.gp1d(:),'double',endian);
fwrite(fileID,master.gw1d(:),'double',endian);

