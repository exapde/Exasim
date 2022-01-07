function [dmd,elem2cpu] = mkpartition1(facecon,f,t2f,t2t,t,nproc,npe,bcm,perm)

ne = size(t,1);
elem2cpu = partition(t,ne,nproc);
for i = 1:nproc
    intelem = find(elem2cpu == (i-1)); % elements in subdomain i
    elem = neighboringelements(t2t, intelem); % all elements connected to elements in subdomain i         
    extelem = setdiff(elem,intelem); % elements outside subdmain i    

    elem = neighboringelements(t2t, extelem); % all elements connected to extelem       
    bndelem = intersect(elem,intelem);  % boundary elements in subdmain i
    outelem = setdiff(elem,[intelem; extelem]); %  elements outside subdmain i            
    
    intelem = setdiff(intelem,bndelem); % interior elements
    dmd{i}.elempart = [intelem; bndelem; extelem; outelem]; % partitions of elements
    dmd{i}.elempartpts = [length(intelem) length(bndelem) length(extelem) length(outelem)];
    nelem = dmd{i}.elempartpts;
    
    recvelem = [extelem; outelem]; % elements received from neighbors
    dmd{i}.elemrecv = [elem2cpu(recvelem)+1 (sum(nelem(1:2))+1:sum(nelem))' recvelem];    
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');  
    dmd{i}.nbsd = unique(dmd{i}.elemrecv(:,1)); % neighboring subdomains
    
    elems = dmd{i}.elempart(1:sum(nelem(1:3))); % [intelem; bndelem; extelem]
    faces = t2f(elems,:); % all faces connected to elems
    faces = unique(faces(:)); % remove duplicate faces
    elemf = f(faces,end-1:end); % elements connected to faces
    bndf = find(elemf(:,2)<=0); % boundary faces    
    intf = setdiff((1:length(faces))',bndf); % interior and exterior faces   
    
    intfaces = t2f(intelem,:); % all faces connected to interior elements
    intfaces = unique(intfaces(:)); % remove duplicate faces
    intfaces = setdiff(intfaces,faces(bndf));
    bndfaces = t2f([bndelem; extelem],:); % all faces connected to boundary elements
    bndfaces = unique(bndfaces(:)); % remove duplicate faces
    %bndfaces = setdiff(bndfaces,faces(bndf));
    insfaces = setdiff(intfaces,bndfaces); % faces within interior elements
%     cmnfaces = setdiff(intfaces,insfaces); % faces between interior elements and boundary elements
%     rmnfaces = setdiff(faces(intf),intfaces); 
%     dmd{i}.facepart = [insfaces; cmnfaces; rmnfaces]; 
%     dmd{i}.facepartpts = [length(insfaces) length(cmnfaces) length(rmnfaces)]; % 
%     dmd{i}.facepartbnd = [0 0 0];
%     dmd{i}.nfbind = cumsum(dmd{i}.facepartpts);        
    rmnfaces = setdiff(faces(intf),insfaces); 
    dmd{i}.facepart = [insfaces; rmnfaces]; 
    dmd{i}.facepartpts = [length(insfaces) length(rmnfaces)]; % 
    dmd{i}.facepartbnd = [0 0];
    dmd{i}.nfbind = cumsum(dmd{i}.facepartpts);        
            
    a = -unique(elemf(bndf,2)); % boundary indices
    bcn = unique(bcm(a));
    nbc = length(bcn); 
    for j=1:nbc % for each boundary condition bcn(j)
        bj = find(bcm==bcn(j)); % find all boundaries that have condition bcn(j)
        bk = [];
        % find all faces that have boundary condition bcn(j)
        for k = 1:length(bj)
            bk = [bk; find(elemf(bndf,2)==-bj(k))];                        
        end        
        % add them to facepart 
        dmd{i}.facepart = [dmd{i}.facepart; faces(bndf(bk))];
        dmd{i}.facepartpts = [dmd{i}.facepartpts length(bk)];
        dmd{i}.facepartbnd = [dmd{i}.facepartbnd bcn(j)];
    end
    if (sum(dmd{i}.facepartpts)-sum(dmd{i}.facepartpts(1:2)))~=length(bndf)
        error('something wrong');
    end    
    
    % break faces into groups so that all faces in one group do not share the same dofs
    mf = cumsum([0 dmd{i}.facepartpts]);  
    bcn = dmd{i}.facepartbnd;
    fcolor = [];
    newmf = 0;
    newbcn = [];
    for j = 1:length(mf)-1 % for each block of faces
        % break the block into smaller blocks        
        [acolor,ncolor] = colorface(f,t2f,perm,dmd{i}.facepart((mf(j)+1):1:mf(j+1))');
        fcolor = [fcolor acolor];
        newmf = [newmf ncolor(2:end)+mf(j)];        
        newbcn = [newbcn bcn(j)*ones(1,length(ncolor(2:end)))];        
    end
    nf = length(faces);
    if (length(unique(fcolor))~=nf) || (newmf(end)~=nf)
        error('something wrong');
    end        
    dmd{i}.facepart = fcolor;    
    dmd{i}.facepartpts = newmf(2:end)-newmf(1:end-1);
    dmd{i}.facepartbnd = newbcn;        
    
    % compute local facecon
    ninf = sum(dmd{i}.facepartpts); % number of faces 
    elemf = f(dmd{i}.facepart(1:ninf),end-1:end);    
    dmd{i}.facecon = facecon(:,:,dmd{i}.facepart(1:ninf));    
    for j = 1:ninf % for each face j       
        fj = dmd{i}.facepart(j); % global face
        ej = elemf(j,end-1:end); % two global elements
        l1 = find(dmd{i}.elempart==ej(1)); % local element for ej(1)
        dmd{i}.facecon(:,1,j) = dmd{i}.facecon(:,1,j) - (ej(1)-l1)*npe;
        if ej(2)>0 % interior face
            l2 = find(dmd{i}.elempart==ej(2)); % local element for ej(2)
            dmd{i}.facecon(:,2,j) = dmd{i}.facecon(:,2,j) - (ej(2)-l2)*npe;
        else % boundary face
            dmd{i}.facecon(:,2,j) = dmd{i}.facecon(:,1,j);
        end
    end    
        
    % local elemen-to-face mapping
    maxf = max(dmd{i}.facepart);
    fmap = zeros(1,maxf);    
    fmap(dmd{i}.facepart) = 1:length(dmd{i}.facepart); % global-to-local face mapping    
    dmd{i}.t2f = fmap(t2f(elems,:)); % get local face indices     
    
    % local face-to-element mapping
    maxe = max(dmd{i}.elempart);
    emap = zeros(1,maxe);    
    emap(dmd{i}.elempart) = 1:length(dmd{i}.elempart); % global-to-local element mapping
    dmd{i}.f = f(dmd{i}.facepart,:);
    dmd{i}.f(:,end-1) = emap(f(dmd{i}.facepart,end-1));    % get local element indices
    ind = f(dmd{i}.facepart,end)>0;  % exclude boundary faces    
    dmd{i}.f(ind,end) = emap(f(dmd{i}.facepart(ind),end)); % get local element indices        
end

for i = 1:nproc
    nei = sum(dmd{i}.elempartpts(1:2));
    elempart{i} = dmd{i}.elempart(1:nei);
end
save dmd.mat elempart;

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc    
    dmd{k}.elemsend = [];
end
for i = 1:nproc            
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        [~,tm(:,2)] = ismember(tm(:,3),dmd{k}.elempart);
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end
% % This is required to sort by entries in first column
%  for k = 1:nproc
%     dmd{k}.elemsend = unique(dmd{k}.elemsend,'rows');        
%  end

for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        dmd{i}.elemsendpts(j) = length(find(dmd{i}.elemsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.elemrecvpts(j) = length(find(dmd{i}.elemrecv(:,1)==dmd{i}.nbsd(j)));
    end
    dmd{i}.elemsend = dmd{i}.elemsend(:,2);
    dmd{i}.elemrecv = dmd{i}.elemrecv(:,2);
end

function nbelem = neighboringelements(t2t, elem)

t2te = t2t(elem,:);            
nbelem = unique(t2te(:));
if nbelem(1) == 0
    nbelem(1) = [];
end
      
function [epart, npart] = partition(t2f,ne,np,weightingElements,elcon)

% Expected formats:
% t2f: ne x nve [2D array]
% elcon: (npf*nfe) x ne [2D array]

if nargin<4; weightingElements = 0; end
if nargin<5; elcon = []; weightingElements = 0; end

if size(t2f,2) == ne; t2f = t2f'; end
if size(elcon,1) == ne; elcon = elcon'; end

% current directory
cdir = pwd;

% move to directory that contains metis programs
if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'exasim') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'EXASIM'))
    up = up+1;
end
if ismac
    cd(strcat(ps(1:is(end-up)),'metis/mac'));
elseif isunix
    cd(strcat(ps(1:is(end-up)),'metis/linux'));
end

% Generate a temporary file to be used in METIS
if weightingElements == 0 || isempty(elcon)
    disp('Writing input files for METIS...');
    dlmwrite('temp.txt', ne, 'delimiter', ' ','precision',10);
    dlmwrite('temp.txt', t2f, '-append', 'delimiter', ' ','precision',10);
else
    disp('Computing element weights for METIS mesh partition...');
    elcon_sorted = sort(elcon,1);
    elcon_sorted = elcon_sorted';
    elcon_sorted(:,end+1) = elcon_sorted(:,end);
    elcon_sorted = elcon_sorted';
    [a,b] = find((elcon_sorted(2:end,:)-elcon_sorted(1:end-1,:)) <= 0);
    [~,i_tmp,~] = unique(b);
    if length(i_tmp) ~= ne; error('Something wrong.'); end
    weights = a(i_tmp).^2;     % We assume cost is O(ncf^2). Actual cost of matrix assembly is O(1) + O(ncf) + O(ncf^2). ncf is the number of unique traced nodes on the element faces
    % Note: We take wieght = ncf^2 to balance the cost of matrix-vector product and
    %       preconditioner solve among processors, as these operations are proportional to the square of
    %       the number of nodes in the 0- and 1-level overlap node subdomain, respectively.
    %       This is so since we make the domain decomposition for nodes such that (no. nodes
    %       in 0-level overlap DD / no. elements in 0-level overlap DD) is approximtely constant for all processors
    % Also, we note that if a few processors do much less work than the
    %       average is much better than if a few processors do much more work
    %       that the average.
    
    disp('Writing input files for METIS...');
    dlmwrite('temp.txt', [ne, 1], 'delimiter', ' ','precision',10);
    dlmwrite('temp.txt', [weights(:), t2f], '-append', 'delimiter', ' ','precision',10);
end

disp('Calling METIS and reading output files...');
if isempty(find(t2f(:,end)==-1, 1))==0  % fix -1
    fin = fopen('temp.txt');
    fout = fopen('output.txt','w');
    while ~feof(fin)
       s = fgetl(fin);
       s = strrep(s, ' -1', '');              
       fprintf(fout,'%s\n',s);       
    end
    fclose(fin);
    fclose(fout);
        
    % call mpmetis
    str = ['!./mpmetis output.txt ' num2str(np)];
    eval(str);

    % get mesh partitioning data
    str = ['output.txt.epart.' num2str(np)];
    epart = textread(str,'%d');

    % get node partitioning data
    str = ['output.txt.npart.' num2str(np)];
    npart = textread(str,'%d');

    % remove files
    delete('temp.txt');
    delete('output.txt');
    str = ['output.txt.epart.' num2str(np)];
    delete(str);
    str = ['output.txt.npart.' num2str(np)];
    delete(str);

    % move back to current directory
    cd(cdir);
else
    % call mpmetis
    str = ['!./mpmetis temp.txt ' num2str(np)];
    eval(str);

    % get mesh partitioning data
    str = ['temp.txt.epart.' num2str(np)];
    epart = textread(str,'%d');

    % get node partitioning data
    str = ['temp.txt.npart.' num2str(np)];
    npart = textread(str,'%d');

    % remove files
    delete('temp.txt');
    str = ['temp.txt.epart.' num2str(np)];
    delete(str);
    str = ['temp.txt.npart.' num2str(np)];
    delete(str);

    % move back to current directory
    cd(cdir);    
end
