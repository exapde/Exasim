function [p,t] = quartercirclemesh(n)

    if mod(n,2) ~= 0
        error('need even number of subdivisions')
    end

    n = n/2;

    [p, t] = squaremesh(n, n, 0, 1);
    
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
    % for ie = 1:size(t,2)
    %    plot(p2(1,[t(:,ie);t(1,ie)],1), p2(2,[t(:,ie);t(1,ie,1)]), 'b-', 'LineWidth', 2); hold on;
    % end
    
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
    
    % for ie = 1:size(t,2)
    %    plot(p3(1,[t(:,ie);t(1,ie)],1), p3(2,[t(:,ie);t(1,ie,1)]), 'k-', 'LineWidth', 2); hold on;
    % end
    
    np = size(p1,2);
    p = [p1, p2, p3];
    t = [t, t+np, t+2*np];

    p = p';
    t = t';
    
    % Remove duplicated nodes:
    snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
    [foo,ix,jx]=unique(round(p/snap)*snap,'rows');
    p=p(ix,:);
    t=jx(t);
    
    % Remove nodes that are not contained in t:
    [pix,ix,jx]=unique(t);
    t=reshape(jx,size(t));
    p=p(pix,:);
   
    p = p';
    t = t';

    lp = 1:size(p,2);
    ind = lp(p(1,:).^2 + p(2,:).^2 > 0.999);
    jnd = lp(p(1,:).^2 + p(2,:).^2 <= 0.999);
    lp = [ind, jnd];

    p = p(:,lp);
    lq = lp*0;
    lq(lp(:)) = 1:size(lp,2);
    t = lq(t);

%     clf;
%     for ie = 1:size(t,2)
%        plot(p(1,[t(:,ie);t(1,ie)],1), p(2,[t(:,ie);t(1,ie,1)]), 'r-', 'LineWidth', 2); hold on;
%     end
%     axis equal;
%     axis tight;
%     pause
end
