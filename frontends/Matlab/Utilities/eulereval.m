function sca = eulereval(u,str,gam,mach)
%EULERVAL Calculates derived quantities for the Euler equation variables.
%   SCA=EULERVAL(U,STR,GAM)
%
%      UP(npl,4,nt):   np plus states
%      STR:            String used to specify requested quantity
%                      - STR: 'r' Density
%                      - STR: 'u' u_x velocity
%                      - STR: 'v' u_y velocity
%                      - STR: 'p' Density
%                      - STR: 'M' Density
%                      - STR; 's' Entropy
%      GAM:            Value of Gamma
%      SCA(npl,4,nt):  Scalar field requested by STR 
%
if strcmp(str,'r')
    sca = u(:,1,:);
elseif strcmp(str,'p')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    sca = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
elseif strcmp(str,'c')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = sqrt(gam*p./u(:,1,:));
elseif strcmp(str,'M')
    u(:,1,:) = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    u2 = sqrt(uv.^2+vv.^2);
    p = abs((gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv)));
    %p = max(p,0.3);
    sca = u2./sqrt(gam*p./u(:,1,:));
elseif strcmp(str,'s')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = p./(u(:,1,:).^gam);
elseif strcmp(str,'u')
    sca = u(:,2,:)./u(:,1,:);
elseif strcmp(str,'v')
    sca = u(:,3,:)./u(:,1,:);
elseif strcmp(str,'c2')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    sca = gam*(gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv))./u(:,1,:);
elseif strcmp(str,'ss')
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    c2  = gam*(gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv))./u(:,1,:);

    r   = u(:,1,:);
    rx  = u(:,5,:);
    rux = u(:,6,:); 
    ry  = u(:,9,:);
    rvy = u(:,11,:);

    r1   = 1./r;

    ux  = (rux - rx.*uv).*r1;
    vy  = (rvy - ry.*vv).*r1;

    sca = (ux+vy)./sqrt(c2);
elseif strcmp(str,'t')
    r  = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p  = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = p./((gam-1)*r);
elseif strcmp(str,'h')
    r  = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p  = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = u(:,4,:)./u(:,1,:)+p./r;    
else
    error('unknonw case');
end