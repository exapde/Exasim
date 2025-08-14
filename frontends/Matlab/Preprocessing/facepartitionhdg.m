function dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc)

%nfe = size(f,1);

for i = 1:nproc
    disp(['face partition ' num2str(i)]); 
        
    % on [intelem bndelem extelem]
    fi = f(:,dmd{i}.elempart); % fix bug here
    dmd{i}.bf = 0*f(:,dmd{i}.elempart);
    for j=1:length(bcm)
      ind = fi==j;
      dmd{i}.bf(ind) = bcm(j);
    end      
    %f2t = mkf2e(t(:,dmd{i}.elempart),elemtype,dim);
    f2t = mkf2t(t(:,dmd{i}.elempart),elemtype,dim);    
    
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
        
%     if nproc>1
%       ne1 = sum(dmd{i}.elempartpts(1:2)); % number of elements in the subdomain
%       f2tin = f2t(:,ina);
%       in1 = find(f2tin(1,:)>ne1);
%       in2 = find(f2tin(3,:)>ne1);
%       in = unique([in1 in2]);      % interdomain faces    
%       ina = [in setdiff(ina, in)]; % [interdomain faces, interior faces]
%       dmd{i}.numinterfacefaces = length(in);
%       t1 = unique(f2tin(1,in2));
%       t2 = unique(dmd{i}.elemsend);
%       if max(abs(t1(:)-t2(:))) > 1e-8        
%         error("domain decomposition is incorrect");
%       end
%     end
    
%     in1
%     in2
%     in
%     inb
%     ina
%     size(ina)
%     size(unique(ina))
%     size(inb)
%     size(f2t)
%     pause
    
    % for i = 1:length(inb)
    %   b = inb(i);
    %   e1 = f2t(1,b);
    %   l1 = f2t(2,b);
    %   fb(i) = fi(l1, e1);
    % end    
    
%     inc = sub2ind(size(fi), f2t(2,inb), f2t(1,inb));
%     fb = fi(inc); % list of boundary indices    
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
        
    bfa = zeros(1, length(inb));      
    for j = 1:length(inb)      
      e1 = f2t(1,inb(j));
      l1 = f2t(2,inb(j));
      bfa(j) = dmd{i}.bf(l1, e1);
    end    
    [sortedbfa, ind] = sort(bfa);
    [facepartbnd, ~, ic] = unique(sortedbfa); %nbc = length(facepartbnd);    
    facepartpts = accumarray(ic, 1); 
    dmd{i}.facepartbnd = [0 facepartbnd];
    dmd{i}.facepartpts = [length(ina) facepartpts(:)'];
        
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

