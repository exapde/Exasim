function meshp = mkmesh_parallel(mesh,epart)
% MESHP = MKMESH_PARALLEL(MESH,EPART) 
% This function generates mesh structure MESHP for each subdomain according to the domain
% partition EPART. MESHP also includes subdomain-to-subdomain connectivities.
%
% MESH  : mesh structure for the whole domain
% EPART : domain partition from metis 
% MESHP : mesh structure for each subdomain

% element-to-element connectivity for the whole domain
t2t = mkt2t(mesh.t,mesh.elemtype); 

% domain partition from metis 
epart = epart+1; % plus one for one base

np = max(epart); % number of subdomains     
npf = size(mesh.perm,1); % number of points on one face
nfe = size(mesh.perm,2); % number of faces on one element

part = cell(np,1);  % ordered element partition for each subdomain
fpart = cell(np,1); % ordered face partition for each subdomain
meshp = cell(np,1); % mesh structure for each subdomain
f = cell(np,1);     % face-to-element connectivity for each subdomain
t2f = cell(np,1);   % element-to-face connectivity for each subdomain

% loop through each subdomain
for i=1:np
    ind = find(epart==i); % list of elements in subdomain i               
    ti = mesh.t(ind,:);   % element-to-vertex connectivity of subdomain i
    t2ti = mkt2t(ti,mesh.elemtype); % element-to-element connectivity of subdomain i
    tj = []; % find boundary elements and interface elements 
    for j=1:size(t2ti,2)
        tj = [tj; find(t2ti(:,j)==0)];
    end
    tb = unique(tj); % remove duplications
    tk = [];  % list of interface elements
    for j=1:length(tb)
        e = ind(tb(j)); % element e
        en = t2t(e,:);  % neighboring elements of element e
        ek = en>0;      % remove zero from en
        em = en(ek);    % neighboring elements of element e
        ei = ismember(em,ind);
        if all(ei)==0 % if any neighboring elements does not belong to subdomain i
            tk = [tk; tb(j)]; % then add this interface element to the list
        end
    end    
    tb = [tk; setdiff(tb,tk)]; % reorder tb to put interface elements first 
    tm = setdiff(1:length(ind),tb); % list of interior elements        
    % ordered element partition
    part{i} = [ind(tb); ind(tm)]; % [interface, boundary, interior]   
    
    % generate face-to-element and element-to-face connectivities for subdomain i 
    [f{i},t2f{i}] = mkt2f(mesh.t(part{i},:),mesh.elemtype);     
        
    % compute face partition
    nf = size(f{i},1);    
    fpart{i} = zeros(nf,1); 
    for n=1:nf % loop through each face
        p = sort(f{i}(n,1:end-2)); % nodes of face n
        % find a face (in mesh.f) that matches face n  
        [a,b] = ismember(p,sort(mesh.f(:,1:end-2),2),'rows');    
        if a==1 % there is always a match
            fpart{i}(n) = b; % put the global face b to the list
        else % if there is no match then it must be a bug
            error('something wrong');
        end
    end
    
    % fixing local f to match with global f
    for n=1:nf
        b = fpart{i}(n); % global face b corresonds to local face n
        f{i}(n,1:end-2) = mesh.f(b,1:end-2); % set nodes of face b to nodes of face n   
        e1 = part{i}(f{i}(n,end-1)); % first global element contains face n           
        e2 = mesh.f(b,end-1); % first global element contains face b
        if (e1~=e2) && f{i}(n,end)>0 
            %if e1 does not match e2 then swap two elements containing face n
            f{i}(n,end-1:end) = f{i}(n,end:-1:end-1);
        end         
    end
    
    % check fpart and t2f
    e = fpart{i}(t2f{i})-mesh.t2f(part{i},:);
    if max(abs(e(:)))>0
        error('something wrong')
    end
    
    % generate mesh structure for subdomain i
    meshp{i}.nd = mesh.nd;
    meshp{i}.elemtype = mesh.elemtype;
    meshp{i}.nodetype = mesh.nodetype;
    meshp{i}.porder = mesh.porder;
    meshp{i}.perm = mesh.perm;
    meshp{i}.tcurved = mesh.tcurved(part{i});
    meshp{i}.p = mesh.p; 
    meshp{i}.t = mesh.t(part{i},:);  % element-to-node connectivity for subdomain i
    meshp{i}.f = [f{i} 0*f{i}(:,1)]; % face-to-element connectivity for subdomain i
    meshp{i}.t2f = t2f{i}; % element-to-face connectivity for subdomain i      
    [meshp{i}.f2f, meshp{i}.fe1, meshp{i}.fe2] = mkf2f(f{i}, t2f{i}); % face-to-face connectivity for subdomain i
    meshp{i}.elcon = elconnectivities(meshp{i}); % element connectivities for subdomain i   
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe length(ind)]);
    meshp{i}.dgnodes = mesh.dgnodes(:,:,part{i}); % dg nodes for subdomain i
    meshp{i}.bf = mesh.bf(:,part{i});             % boundary faces for subdomain i
    meshp{i}.eint = length(tk); % number of interface elements for subdomain i
    meshp{i}.epart = part{i};   % global indices of elements of subdomain i
    meshp{i}.f = f{i};          % face-to-element connectivity for subdomain i
    meshp{i}.fpart = fpart{i};  % global indices of faces of subdomain i
    meshp{i}.plocal = mesh.plocal;
    meshp{i}.tlocal = mesh.tlocal;
    meshp{i}.plocfc = mesh.plocfc;
    meshp{i}.tlocfc = mesh.tlocfc;  
        
    % compute subdomain-to-subdomain connectivity
    for j=1:np
        meshp{i}.emap{j} = []; % element connectivity between subdmomains
        meshp{i}.fmap{j} = []; % face connectivity between subdmomains
    end    
    meshp{i}.pmap = zeros(1,np); % processor connectivity between subdmomains
    bf = find(meshp{i}.f(:,end) == 0); % boundary faces
    for m = 1:length(bf) % for each boundary face
        fm = bf(m); % local boundary face fm
        ei = meshp{i}.f(fm,end-1); % local element ei which contains fm 
        if ei <= meshp{i}.eint % check if em is an interface element                       
            eg = part{i}(ei);  % global element eg which matches ei
            fg = fpart{i}(fm); % global face fg which matches fm            
            if mesh.f(fg,end)>0 % check if fg is an interface boundary
                en = setdiff(mesh.f(fg,end-1:end),eg); % neighboring element en of eg
                j  = epart(en); % subdomain j which contains en
                meshp{i}.pmap(j) = 1;
                ja = find(mesh.t2f(en,:)==fg); % position of face fg in element en                         
                jj = [ja setdiff(1:nfe,ja)]; % [first face fg, then the remaining faces]        
                % [local element ei, positions of the faces in element en]
                meshp{i}.emap{j} = [meshp{i}.emap{j}; fg ei jj]; 
                ia = find(mesh.t2f(eg,:)==fg); % position of face fg in element eg 
                ii = [ia setdiff(1:nfe,ia)]; % [first face fg, then the remaining faces]                      
                % local faces of local element ei
                meshp{i}.fmap{j} = [meshp{i}.fmap{j}; [fg meshp{i}.t2f(ei,ii)]];
            end
        end
    end
    for j=1:np
        if meshp{i}.pmap(j) == 1
            % sort the global faces to obtain the indices
            [~,ind] = sort(meshp{i}.fmap{j}(:,1));
            % make sure to match meshp{i}.emap{j} to meshp{j}.emap{i} 
            meshp{i}.emap{j} = meshp{i}.emap{j}(ind,2:end);
            % make sure to match meshp{i}.fmap{j} to meshp{j}.fmap{i} 
            meshp{i}.fmap{j} = meshp{i}.fmap{j}(ind,2:end);            
        end
    end
    meshp{i}.pmap = find(meshp{i}.pmap==1);
end

 
% pmap = zeros(np,np);
% emap = cell(np,np);
% fmap = cell(np,np);
% gmap = cell(np,np);
% lmap = cell(np,np);
% for i=1:np
%     for j=1:np
%         emap{i}{j} = [];
%         fmap{i}{j} = [];
%         gmap{i}{j} = [];
%         lmap{i}{j} = [];
%     end
% end
% 
% % tpart = zeros(size(mesh.t,1),1);
% % for i=1:np
% %     tpart(part{i}) = i;    
% % end
% % for i=1:np
% %     for j=1:np
% %         meshp{i}.emap2{j} = [];
% %         meshp{i}.fmap2{j} = [];
% %     end
% % end
% % for i = 1:np % for each subdomain i
% %     for j=1:np
% %         meshp{i}.emap2{j} = [];
% %         meshp{i}.fmap2{j} = [];
% %     end    
% %     meshp{i}.pmap2 = zeros(1,np);
% %     bf = find(meshp{i}.f(:,end) == 0); % boundary faces
% %     for m = 1:length(bf) % for each boundary face
% %         fm = bf(m); % local boundary face fm
% %         ei = meshp{i}.f(fm,end-1); % local element ei which contains fm 
% %         if ei <= meshp{i}.eint % check if em is an interface element                       
% %             eg = part{i}(ei);  % global element eg which matches em
% %             fg = fpart{i}(fm); % global face fg which matches fm            
% %             if mesh.f(fg,end)>0 % check if fg is an interface boundary
% %                 en = setdiff(mesh.f(fg,end-1:end),eg); % neighboring element en of eg
% %                 j  = epart(en); % subdomain j which contains en
% %                 meshp{i}.pmap2(j) = 1;
% %                 ja = find(mesh.t2f(en,:)==fg); % position of face fg in element en                         
% %                 jj = [ja setdiff(1:nfe,ja)]; % [first face fg, then the remaining faces]        
% %                 % [local element ei, positions of the faces in element en]
% %                 meshp{i}.emap2{j} = [meshp{i}.emap2{j}; ei jj]; 
% %                 ia = find(mesh.t2f(eg,:)==fg); % position of face fg in element eg 
% %                 ii = [ia setdiff(1:nfe,ia)]; % [first face fg, then the remaining faces]                      
% %                 % local faces of local element ei
% %                 meshp{i}.fmap2{j} = [meshp{i}.fmap2{j}; meshp{i}.t2f(ei,ii)];
% %             end
% %         end
% %     end
% %     meshp{i}.pmap2 = find(meshp{i}.pmap2==1);
% % end
% 
% nf = size(mesh.f,1);
% for n=1:nf % for each global face n
%     fn = mesh.f(n,end-1:end); % two elements sharing face n
%     if fn(2)>0 % if face n is interior, then do
%         % find subdomain i that contains element fn(1)
%         for i=1:np % for each subdomain
%             ei = find(part{i}==fn(1)); % find local element ei in subdomain i to match element fn(1)
%             if isempty(ei)==0
%                 break;
%             end
%         end
%         % find subdomain j that contains element fn(2)
%         for j=1:np % for each subdomain           
%             ej = find(part{j}==fn(2)); % find local element ej in subdomain j to match element fn(2)
%             if isempty(ej)==0
%                 break;
%             end
%         end        
%         if i~=j % if subdomain i differs from subdomain j, then it is an interface face
%             pmap(i,j) = 1; % connect subdomain i to subdomain j
%             pmap(j,i) = 1; % connect subdomain j to subdomain i
%             mi = find(fpart{i}==n); % find a local face mi in subdomain i to match global face n
%             mj = find(fpart{j}==n); % find a local face mj in subdomain j to match global face n
%             ia = find(t2f{i}(ei,:)==mi); % position of face mi in element ei 
%             ii = [ia setdiff(1:nfe,ia)]; % [first face mi, then the remaining faces]                      
%             ja = find(t2f{j}(ej,:)==mj); % position of face mj in element ej                         
%             jj = [ja setdiff(1:nfe,ja)]; % [first face mj, then the remaining faces]        
%             emap{i}{j} = [emap{i}{j}; [ei jj]]; % [local element ei, positions of the faces in element ej]
%             emap{j}{i} = [emap{j}{i}; [ej ii]]; % [local element ej, positions of the faces in element ei]
%             fmap{i}{j} = [fmap{i}{j}; t2f{i}(ei,ii)]; % local faces of local element ei
%             fmap{j}{i} = [fmap{j}{i}; t2f{j}(ej,jj)]; % local faces of local element ej                         
%             gmap{i}{j} = [gmap{i}{j}; [fn(1) fn(2)]];
%             gmap{j}{i} = [gmap{j}{i}; [fn(2) fn(1)]];            
%             lmap{i}{j} = [lmap{i}{j}; [ei ej]];
%             lmap{j}{i} = [lmap{j}{i}; [ej ei]];            
%         end
%     end
% end
% 
% % for i=1:np
% %     a1(i) = max(abs(meshp{i}.pmap -find(pmap(i,:)==1)));    
% %     for j=1:np
% %         if isempty(meshp{i}.emap{j})==0
% %             e = max(abs(meshp{i}.emap{j}(:)-emap{i}{j}(:)));
% %             if e>0
% %                 meshp{i}.emap{j}
% %                 emap{i}{j}
% %                 pause
% %             end
% %             e = max(abs(meshp{i}.fmap{j}(:)-fmap{i}{j}(:)));
% %             if e>0
% %                 meshp{i}.fmap{j}
% %                 fmap{i}{j}
% %                 pause
% %             end            
% %         end
% %     end
% % end
% 
% 
% % subdomain-to-subdomain connectivities
% % for i=1:np    
% %     meshp{i}.fint = 0; % number of faces on the interface for each subdomain   
% %     for j=1:np
% %         meshp{i}.pmap = find(pmap(i,:)==1); % processor connectivity between subdmomains
% %         meshp{i}.emap{j} = emap{i}{j};      % element connectivity between subdmomains
% %         meshp{i}.fmap{j} = fmap{i}{j};      % face connectivity between subdmomains
% %         if pmap(i,j)==1
% %             meshp{i}.fint = meshp{i}.fint + size(fmap{i}{j},1);
% %         end
% %     end
% % end
% 
% if mesh.nd==2
%     % plot the mesh structures 
%     bcol = [1 1 0;... % yellow 
%             1 0 1;... % magenta
%             0 1 1;... % cyan
%             1 0 0;... % red
%             0 1 0;... % green 
%             0 0 1;... % blue
%             1,0.4,0.6;...
%             0.4,0.6,1;...
%            ];
%    
%     figure(1); clf;
%     hold on;            
%     for i=1:np
%         simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
%     end
%     for i=1:np-1
%         for j=i+1:np
%             e = gmap{i}{j}(:);
%             for it=1:length(e)
%                 pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
%                 txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%                 text(pmid(1),pmid(2),num2str(e(it)),txtpars{:});
%             end        
%         end
%     end
%     hold off;
%     axis equal;      
%     axis tight;
%     axis off;  
% 
%     figure(2); clf;
%     hold on;            
%     for i=1:np
%         simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
%     end
%     for i=1:np-1
%         for j=i+1:np
%             e = gmap{i}{j}(:);
%             l = lmap{i}{j}(:);
%             for it=1:length(e)
%                 pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
%                 txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%                 text(pmid(1),pmid(2),num2str(l(it)),txtpars{:});
%             end        
%         end
%     end
%     hold off;
%     axis equal;      
%     axis tight;
%     axis off;  
%     
%     figure(3); clf;
%     hold on;            
%     for i=1:np
%         simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
%     end
%     for i=1:np
%         for j=1:np
%             if isempty(gmap{i}{j})==0
%                 e = gmap{i}{j}(:,1);
%                 f = fmap{i}{j}(:,1);
%                 for it=1:length(e)
%                     pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
%                     txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%                     %text(pmid(1),pmid(2),num2str(f(it)),txtpars{:});
%                     text(pmid(1)*(0.95+(1.05-0.95)*rand),pmid(2)*(0.95+(1.05-0.95)*rand),num2str(f(it)),txtpars{:});
%                 end    
%             end
%         end
%     end
%     hold off;
%     axis equal;      
%     axis tight;
%     axis off;
%     
%     figure(4); clf;
%     hold on;            
%     for i=1:np
%         simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
%     end
%     for i=1:np
%         for j=1:np
%             if isempty(gmap{i}{j})==0
%                 e = gmap{i}{j}(:,1);
%                 f = fmap{i}{j}(:,1);
%                 g = fpart{i}(f);
%                 for it=1:length(e)
%                     pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
%                     txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%                     text(pmid(1)*(0.95+(1.05-0.95)*rand),pmid(2)*(0.95+(1.05-0.95)*rand),num2str(g(it)),txtpars{:});
%                 end    
%             end
%         end
%     end
%     hold off;
%     axis equal;      
%     axis tight;
%     axis off;
% end
% 
