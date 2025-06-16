function dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc)

nfe = size(f,1);

for i = 1:nproc
    disp(['face partition ' num2str(i)]); 
        
    % on [intelem bndelem extelem]
    fi = f(:,dmd{i}.elempart); % fix bug here
    dmd{i}.bf = 0*f(:,dmd{i}.elempart);
    for j=1:length(bcm)
      ind = fi==j;
      dmd{i}.bf(ind) = bcm(j);
    end      
    f2t = mkf2e(t(:,dmd{i}.elempart),elemtype,dim);
        
    % only on [intelem bndelem]
    nelem = dmd{i}.elempartpts;
    if numel(nelem)>=2
        ne3 = sum(nelem(1:2));
        ind = f2t(1,:)<=ne3;
        f2t = f2t(:,ind);
    end
            
    % reorder so that boundary faces are last
    ina = find(f2t(3,:)>0); % interior faces
    inb = find(f2t(3,:)==0); % boundary faces
    inc = sub2ind(size(fi), f2t(2,inb), f2t(1,inb));
    fb = fi(inc); % list of boundary indices
          
    if nproc>1
      ne1 = sum(dmd{i}.elempartpts(1:2)); % number of elements in the subdomain
      f2tin = f2t(:,ina);
      in1 = find(f2tin(1,:)>ne1);
      in2 = find(f2tin(3,:)>ne1);
      in = unique([in1 in2]);      % interface faces    
      ina = [in setdiff(ina, in)]; % [interface faces, interior faces]
      dmd{i}.numinterfacefaces = length(in);
      t1 = unique(f2tin(1,in2));
      t2 = unique(dmd{i}.elemsend);
      if max(abs(t1(:)-t2(:))) > 1e-8        
        error("domain decomposition is incorrect");
      end
    end
    
    fa = unique(fb); % boundary indices    
    bcn = unique(bcm(fa)); % a list of boundary conditions
    nbc = length(bcn);

    dmd{i}.facepartpts = length(ina);
    dmd{i}.facepartbnd = 0;
    ind = zeros(1,length(fb));
    m = 1;
    for j=1:nbc % for each boundary condition bcn(j)
        bj = find(bcm==bcn(j)); % find all boundaries that have condition bcn(j)
        n = 0;
        for k = 1:length(bj) % for each boundary that has condition bcn(j)
            ii = find(fb == bj(k)); % indices of the boundary bj(k)
            l = length(ii);
            n = n + l;
            ind(m:(m+l-1)) = ii;
            m = m + l;
        end
        dmd{i}.facepartpts = [dmd{i}.facepartpts n];
        dmd{i}.facepartbnd = [dmd{i}.facepartbnd bcn(j)];
    end    
    
    % [interface faces, interior faces, boundary faces]
    f2t = f2t(:,[ina inb(ind)]);
    dmd{i}.f2t = f2t;    
    
    % fix f2t to ensure DOF consistency across interior faces
    % the global ID of the FIRST element must be smaller than that of the SECOND element        
    N = size(dmd{i}.f2t,2);
    for face = 1:N % loop over each face
      e1 = dmd{i}.f2t(1,face);
      e2 = dmd{i}.f2t(3,face);
      if e2>0 % if the face is an interior face
        g1 = dmd{i}.elempart(e1);
        g2 = dmd{i}.elempart(e2);
        if (g2<g1) 
          % ensure g2 > g1      
          tm = dmd{i}.f2t(3:4,face);
          dmd{i}.f2t(3:4,face) = dmd{i}.f2t(1:2,face);
          dmd{i}.f2t(1:2,face) = tm;            
        end
        e1 = dmd{i}.f2t(1,face);
        e2 = dmd{i}.f2t(3,face);
        g1 = dmd{i}.elempart(e1);
        g2 = dmd{i}.elempart(e2);
        if (g2<=g1)
          error("something wrong");
        end
      end
    end     
%     dmd{i}.t2f = mke2f(dmd{i}.f2t);    
%     N = length(dmd{i}.elemsend);
%     for n = 1:N % loop over each interface element/face
%       esend = dmd{i}.elemsend(n);
%       erecv = dmd{i}.elemrecv(n);
%       for j = 1:nfe
%         fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
%         for k = 1:nfe
%           frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
%           if (fsend == frecv)          % if esend and erecv share the same face
%             dmd{i}.facesend(n) = j;    % local ID of the face sent to neighbor
%             dmd{i}.facerecv(n) = k;    % local ID of the face received from neighbor                      
%           end
%         end
%       end                    
%     end         
    
    [dmd{i}.facecon,dmd{i}.elemcon] = faceconnectivity2(t(:,dmd{i}.elempart),dmd{i}.f2t,dim,elemtype,porder);      
end


% nfe = size(f,1);
% 
% % % remove faces in f2t to avoid DOF duplication across interface faces 
% % for i=1:nproc
% %   gsend = dmd{i}.elempart(dmd{i}.elemsend); % global ID of elements sent to neighbors
% %   grecv = dmd{i}.elempart(dmd{i}.elemrecv); % global ID of elements received from neighbors   
% %   dmd{i}.gsendrecv = [gsend(:) grecv(:) 0*grecv(:)];
% % end
% % % 1 -> kept, 0 -> removed
% % dmd{1}.gsendrecv(1:2:end,3) = 1; % kept
% % dmd{1}.gsendrecv(2:2:end,3) = 0; % removed
% % for i=2:length(dmd)
% %   N = length(dmd{i}.elemsend);
% %   token = 1;
% %   for j = 1:N
% %     nb = dmd{i}.nbsend(j); % neighbor ID
% %     if nb < i      
% %       for k = 1:length(dmd{nb}.elemsend)
% %         if (dmd{i}.gsendrecv(j,1) == dmd{nb}.gsendrecv(k,2)) && (dmd{i}.gsendrecv(j,2) == dmd{nb}.gsendrecv(k,1))
% %           dmd{i}.gsendrecv(j,3) = 1 - dmd{nb}.gsendrecv(k,3);        
% %         end
% %       end
% %     else
% %       dmd{i}.gsendrecv(j,3) = token;
% %       token = 1-token;
% %     end    
% %   end
% % end
% 
% for i = 1:nproc
%     disp(['face partition ' num2str(i)]); 
%     
%     % on [intelem bndelem extelem]
%     dmd{i}.bf = f(:,dmd{i}.elempart);
%     for j=1:length(bcm)
%       ind = dmd{i}.bf==j;
%       dmd{i}.bf(ind) = bcm(j);
%     end
%     fi = f(:,dmd{i}.elempart);
%     f2t = mkf2e(t(:,dmd{i}.elempart),elemtype,dim);
% 
%     % only on [intelem bndelem]
%     nelem = dmd{i}.elempartpts;
%     if numel(nelem)>=2
%         ne3 = sum(nelem(1:2));
%         ind = f2t(1,:)<=ne3;
%         f2t = f2t(:,ind);
%     end
% 
% %     dmd{i}.t2f = mke2f(f2t);
% %         
% %     % remove faces in f2t to avoid DOF duplication across interface faces 
% %     N = size(dmd{i}.gsendrecv,1);
% %     dmd{i}.removed = zeros(N,1);  
% %     for n = 1:N % loop over each interface element/face
% %       if (dmd{i}.gsendrecv(n,3)==0) % face to be removed
% %         esend = dmd{i}.elemsend(n); % element sent to neighbor
% %         erecv = dmd{i}.elemrecv(n); % element received from neighbor
% %         for j = 1:nfe
% %           fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
% %           for k = 1:nfe
% %             frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
% %             if (fsend == frecv)          % if esend and erecv share the same face          
% %               dmd{i}.removed(n) = fsend;
% %             end
% %           end
% %         end
% %       end
% %     end
% %     f2t(:,dmd{i}.removed(dmd{i}.removed > 0)) = [];    
%     
%     % reorder so that boundary faces are last
%     ina = find(f2t(3,:)>0); % interior faces
%     inb = find(f2t(3,:)==0); % boundary faces
%     inc = sub2ind(size(fi), f2t(2,inb), f2t(1,inb));
%     fb = fi(inc); % list of boundary indices
%         
%     fa = unique(fb); % boundary indices    
%     bcn = unique(bcm(fa)); % a list of boundary conditions
%     nbc = length(bcn);
% 
%     dmd{i}.facepartpts = length(ina);
%     dmd{i}.facepartbnd = 0;
%     ind = zeros(1,length(fb));
%     m = 1;
%     for j=1:nbc % for each boundary condition bcn(j)
%         bj = find(bcm==bcn(j)); % find all boundaries that have condition bcn(j)
%         n = 0;
%         for k = 1:length(bj) % for each boundary that has condition bcn(j)
%             ii = find(fb == bj(k)); % indices of the boundary bj(k)
%             l = length(ii);
%             n = n + l;
%             ind(m:(m+l-1)) = ii;
%             m = m + l;
%         end
%         dmd{i}.facepartpts = [dmd{i}.facepartpts n];
%         dmd{i}.facepartbnd = [dmd{i}.facepartbnd bcn(j)];
%     end
%     % [interior faces, boundary faces]
%     f2t = f2t(:,[ina inb(ind)]);
%     dmd{i}.f2t = f2t;    
% 
%     % fix f2t to ensure DOF consistency across interior faces
%     % the global ID of the FIRST face must be smaller than that of the SECOND face        
%     N = size(dmd{i}.f2t,2);
%     for face = 1:N % loop over each face
%       e1 = dmd{i}.f2t(1,face);
%       e2 = dmd{i}.f2t(3,face);
%       if e2>0 % if the face is an interior face
%         g1 = dmd{i}.elempart(e1);
%         g2 = dmd{i}.elempart(e2);
%         if (g2<g1) 
%           % ensure g2 > g1      
%           tm = dmd{i}.f2t(3:4,face);
%           dmd{i}.f2t(3:4,face) = dmd{i}.f2t(1:2,face);
%           dmd{i}.f2t(1:2,face) = tm;            
%         end
%         e1 = dmd{i}.f2t(1,face);
%         e2 = dmd{i}.f2t(3,face);
%         g1 = dmd{i}.elempart(e1);
%         g2 = dmd{i}.elempart(e2);
%         if (g2<=g1)
%           error("something wrong");
%         end
%       end
%     end     
%     
%     dmd{i}.t2f = mke2f(dmd{i}.f2t);    
%     N = length(dmd{i}.elemsend);
%     for n = 1:N % loop over each interface element/face
%       esend = dmd{i}.elemsend(n);
%       erecv = dmd{i}.elemrecv(n);
%       for j = 1:nfe
%         fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
%         for k = 1:nfe
%           frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
%           if (fsend == frecv)          % if esend and erecv share the same face
%             dmd{i}.facesend(n) = j;    % local ID of the face sent to neighbor
%             dmd{i}.facerecv(n) = k;    % local ID of the face received from neighbor                      
%           end
%         end
%       end                    
%     end         
%     
%     [dmd{i}.facecon,dmd{i}.elemcon] = faceconnectivity2(t(:,dmd{i}.elempart),dmd{i}.f2t,dim,elemtype,porder);      
% end
% 
% % % remove faces in f2t to avoid DOF duplication across interface faces 
% % for i=1:length(dmd)
% %   N = length(dmd{i}.elemsend);
% %   dmd{i}.removed = zeros(N,1);  
% %   for n = 1:N % loop over each interface element/face
% %     if (dmd{i}.gsendrecv(n,3)==0) % face to be removed
% %       esend = dmd{i}.elemsend(n); % element sent to neighbor
% %       erecv = dmd{i}.elemrecv(n); % element received from neighbor
% %       for j = 1:nfe
% %         fsend = dmd{i}.t2f(j,esend);   % face sent to neighbor
% %         for k = 1:nfe
% %           frecv = dmd{i}.t2f(k,erecv); % face received from neighbor
% %           if (fsend == frecv)          % if esend and erecv share the same face          
% %             dmd{i}.removed(n) = fsend;
% %           end
% %         end
% %       end
% %     end
% %   end
% %   dmd{i}.f2t(:,dmd{i}.removed(dmd{i}.removed > 0)) = [];
% % end
% 
% % for i=1:length(dmd)
% %   [dmd{i}.facecon,dmd{i}.elemcon] = faceconnectivity2(t(:,dmd{i}.elempart),dmd{i}.f2t,dim,elemtype,porder);      
% % end
% 
