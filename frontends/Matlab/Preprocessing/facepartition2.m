function dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc)


for i = 1:nproc
    disp(['face partition ' num2str(i)]); 
    dmd{i}.bf = f(:,dmd{i}.elempart);
    for j=1:length(bcm)
      ind = dmd{i}.bf==j;
      dmd{i}.bf(ind) = bcm(j);
    end
    
    fi = f(:,dmd{i}.elempart);
    f2t = mkf2e(t(:,dmd{i}.elempart),elemtype,dim);

    % only on [intelem bndelem extelem]
    nelem = dmd{i}.elempartpts;
    if numel(nelem)>=3
        ne3 = sum(nelem(1:3));
        ind = f2t(1,:)<=ne3;
        f2t = f2t(:,ind);
    end

    % reorder so that boundary faces are last
    ina = find(f2t(3,:)>0); % interior faces
    inb = find(f2t(3,:)==0); % boundary faces
    inc = sub2ind(size(fi), f2t(2,inb), f2t(1,inb));
    fb = fi(inc); % list of boundary indices
        
    % for i = 1:length(inb)
    %   b = inb(i);
    %   e1 = f2t(1,b);
    %   l1 = f2t(2,b);
    %   fb(i) = fi(l1, e1);
    % end
    % bcn = bcm(fb);
    % [b,ind]=sort(bcn);
    % f2t = f2t(:,[ina inb(ind)]);
    
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
    % [interior faces, boundary faces]
    f2t = f2t(:,[ina inb(ind)]);
    dmd{i}.f2t = f2t;
    dmd{i}.t2f = mke2f(f2t);
    [dmd{i}.facecon,dmd{i}.elemcon] = faceconnectivity2(t(:,dmd{i}.elempart),f2t,dim,elemtype,porder);    
end



