function dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc)


for i = 1:nproc
    disp(['face partition ' num2str(i)]); 
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
    dmd{i}.facecon = faceconnectivity2(t(:,dmd{i}.elempart),f2t,dim,elemtype,porder);    
end


