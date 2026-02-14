function mesh = mkmesh_quartercircle(porder,n)

    if mod(n,2) ~= 0
        error('need even number of subdivisions')
    end

    n = n/2;
    mesh = mkmesh_square(n,n,porder,1,1,1,1,1);
    [mesh.p, mesh.t, mesh.dgnodes] = square2quartercircle(mesh.p, mesh.t, mesh.dgnodes, n);
    mesh.xpe = mesh.plocal;
    mesh.telem = mesh.tlocal;

    figure(1); clf; meshplot(mesh,1);
end

function [p, t, xdg] = square2quartercircle(p, t, xdg, n)
    
    np = size(p,2);    
    t = [t, t+np, t+2*np];

    p = mapping(p, n);
    xdg = mapping2(xdg, n);
    
    % [npe, nd, ne] = size(xdg);
    % xdg = reshape(permute(xdg, [2 1 3]), [nd npe*ne]);
    % xdg = mapping(xdg, n);

    p = p';
    t = t';
    
    % Remove duplicated nodes:
    snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
    [~,ix,jx]=unique(round(p/snap)*snap,'rows');
    p=p(ix,:);
    t=jx(t);
    
    % Remove nodes that are not contained in t:
    [pix,~,jx]=unique(t);
    t=reshape(jx,size(t));
    p=p(pix,:);
   
    p = p';
    t = t';

    % lp = 1:size(p,2);
    % ind = lp(p(1,:).^2 + p(2,:).^2 > 0.999);
    % jnd = lp(p(1,:).^2 + p(2,:).^2 <= 0.999);
    % lp = [ind, jnd];
    % 
    % p = p(:,lp);
    % lq = lp*0;
    % lq(lp(:)) = 1:size(lp,2);
    % t = lq(t);
end

function p = mapping(p, n)

    th = p(1,:)*pi/4.0;
    rmin = 2.*0.5*(0.0 + 0.5./cos(th));
    rmin = 0.25 + 0.5*rmin;
    r =  rmin + (1 - rmin).*p(2,:);
    
    p1 = 0*p;
    p1(1,:) = r.*sin(th);
    p1(2,:) = r.*cos(th);
    
    th = pi/4.0 + p(1,:)*pi/4.0;
    rmin = 0.5./cos(pi/2.0 - th);
    rmin = 0.25 + 0.5*rmin;
    r =  rmin + (1 - rmin).*p(2,:);
    
    p2 = 0*p;
    p2(1,:) = r.*sin(th);
    p2(2,:) = r.*cos(th);
    
    p3 = 0*p;
    
    up = p1(:,1:n+1);
    bot = p(:,1:n+1)*0.5;
    ver = p(2,1:n+1:end);
    
    lef = p(:,1:n+1:end)*0.5;
    rig = flip(p2(:,1:n+1),2);
    hor = p(1,1:n+1);
    
    x1 = (1-hor)'*lef(1,:) + hor'*rig(1,:);
    y1 = (1-hor)'*lef(2,:) + hor'*rig(2,:);
    
    x2 = bot(1,:)'*(1-ver) + up(1,:)'*ver;
    y2 = bot(2,:)'*(1-ver) + up(2,:)'*ver;
    
    xi = - (1-hor)'*lef(1,1)*(1-ver) - (1-hor)'*lef(1,end)*ver - hor'*rig(1,1)*(1-ver) - hor'*rig(1,end)*ver;
    yi = - (1-hor)'*lef(2,1)*(1-ver) - (1-hor)'*lef(2,end)*ver - hor'*rig(2,1)*(1-ver) - hor'*rig(2,end)*ver;    

    p3(1,:) = x1(:) + x2(:) + xi(:);
    p3(2,:) = y1(:) + y2(:) + yi(:);
            
    p = [p1, p2, p3];    
end

function p = mapping2(p, n)

    [npe, nd, ne] = size(p);
    k = sqrt(npe);

    th = p(:,1,:)*pi/4.0;
    rmin = 2.*0.5*(0.0 + 0.5./cos(th));
    rmin = 0.25 + 0.5*rmin;
    r =  rmin + (1 - rmin).*p(:,2,:);
    

    p1 = 0*p;
    p1(:,1,:) = r.*sin(th);
    p1(:,2,:) = r.*cos(th);
    
    th = pi/4.0 + p(:,1,:)*pi/4.0;
    rmin = 0.5./cos(pi/2.0 - th);
    rmin = 0.25 + 0.5*rmin;
    r =  rmin + (1 - rmin).*p(:,2,:);
    
    p2 = 0*p;
    p2(:,1,:) = r.*sin(th);
    p2(:,2,:) = r.*cos(th);
    
    p3 = 0*reshape(permute(p, [2 1 3]), [nd npe*ne]);
    
    up = p1(1:k,:,1:n);
    bot = p(1:k,:,1:n)*0.5;
    %ver = p(1:k,2,1:n:end);
    
    lef = p(1:k:npe,:,1:n:end)*0.5;
    rig = p2(1:k,:,1:n);
    hor = p(1:k,1,1:n);
    
    up = reshape(permute(up, [2 1 3]), [nd k*n]);
    bot = reshape(permute(bot, [2 1 3]), [nd k*n]);
    %ver = reshape(permute(ver, [2 1 3]), [1 k*n]);

    lef = reshape(permute(lef, [2 1 3]), [nd k*n]);
    rig = flip(reshape(permute(rig, [2 1 3]), [nd k*n]), 2);
    hor = reshape(permute(hor, [2 1 3]), [1 k*n]);
    ver = hor;

    x1 = (1-hor)'*lef(1,:) + hor'*rig(1,:);
    y1 = (1-hor)'*lef(2,:) + hor'*rig(2,:);
    
    x2 = bot(1,:)'*(1-ver) + up(1,:)'*ver;
    y2 = bot(2,:)'*(1-ver) + up(2,:)'*ver;
    
    xi = - (1-hor)'*lef(1,1)*(1-ver) - (1-hor)'*lef(1,end)*ver - hor'*rig(1,1)*(1-ver) - hor'*rig(1,end)*ver;
    yi = - (1-hor)'*lef(2,1)*(1-ver) - (1-hor)'*lef(2,end)*ver - hor'*rig(2,1)*(1-ver) - hor'*rig(2,end)*ver;    
    
    p3(1,:) = x1(:) + x2(:) + xi(:);
    p3(2,:) = y1(:) + y2(:) + yi(:);
     
    p3 = permute(reshape(p3, [nd k n k n]), [1 2 4 3 5]);
    p3 = permute(reshape(p3, [nd npe ne]), [2 1 3]);

    p = cat(3, p1, p2);
    p = cat(3, p, p3);        

    % figure(2); clf; hold on;
    % plot(up(1,:), up(2,:), 'o');
    % plot(bot(1,:), bot(2,:), 'o');
    % plot(lef(1,:), lef(2,:), 'o');
    % plot(rig(1,:), rig(2,:), 'o');   
end
