function [fh,fh_udg,fh_uh] = fbou(ib,uinf,nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):             Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

nc = size(udg,2);
[ng,nch] = size(uh);
nd = nch-2;

if nch==5    
    [fh,fh_udg,fh_uh] = fbou3d(ib,uinf,nl,p,udg,uh,param,time);
    return;
end

gam = param{1};
gam1 = gam-1.0;
Re = param{3};
Minf = param{5};
M2   = Minf^2;


u = udg(:,1:nch);
q = reshape(udg(:,nch+1:3*nch),ng,nch,2);

switch ib
    case 1  % Far field
%         [an,anm] = getan(nl,uh,param,0);
%         [An,Anm] = getan(nl,uh,param,1);                
%         fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
% 
%         fh_u = an+An;
%         fh_q = zeros(ng,nch,nch,2);
%         fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
%         
%         fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));      
        
        [An,Anm] = getan2d(nl,uh,param);                
        fh = 0.5*(u + uinf) + 0.5*permute(mapContractK(An,u-uinf,2,3,1,2,[],1),[2 1]) - uh;
        fh_u = 0.5*An;
        for i = 1:nch
          fh_u(:,i,i) = fh_u(:,i,i) + 0.5;
        end
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = 0.5*permute(mapContractK(Anm,u-uinf,[2 4],3,1,2,[],1),[3 1 2]);
        for i = 1:nch
          fh_uh(:,i,i) = fh_uh(:,i,i) - 1.0;
        end
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));   
        
%         fb = tmp(nl,udg,uh,uinf,gam);        
%         fb_u = tmp2(nl,udg,uh,uinf,gam);    
%         fb_u = reshape(fb_u,ng,nch,nch);   
%         fb_uh = tmp3(nl,udg,uh,uinf,gam);     
%         fb_uh = reshape(fb_uh,ng,nch,nch);
%         
%         fh = reshape(fb,[ng nch]);
%         fh_uh = fb_uh;  
%         fh_udg = cat(3,fb_u,zeros(ng,nch,2*nch));         
        
%         max(abs(fh(:)-fb(:)))                        
%         max(abs(fh_u(:)-fb_u(:)))                
%         max(abs(fh_uh(:)-fb_uh(:)))        
    case 2  % Adiabatic Wall     
        uinf = u;
        uinf(:,2:3) = 0;
        
        [fh4,fh4_udg,fh4_uh] = fhat(nl,p,udg,uh,param,time);
        
        fh = uinf - uh;
        fh(:,4) = fh4(:,4);
        
        fh_u = zeros(ng,nch,nch);
        fh_u(:,1,1) = ones(ng,1);        
        fh_u(:,4,:) = fh4_udg(:,4,1:nch);
   
        fh_q = zeros(ng,nch,nch*2);
        fh_q(:,4,:) = fh4_udg(:,4,nch+1:3*nch);
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,:) = fh4_uh(:,4,:);
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 3 % Isothermal wall
        %[fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        %fh_u = fh_udg(:,:,1:nch);
        %fh_q = fh_udg(:,:,nch+1:3*nch);
        
        %one = ones(ng,1);        
        zero = zeros(ng,1);
        fh(:,1)   = u(:,1)-uh(:,1);        
        fh(:,2:3) = -uh(:,2:3);
        fh(:,4) = gam*gam1*M2*uh(:,4)./uh(:,1) - uinf(:,end);
        
        fh_u = zeros(ng,nch,nch);
        fh_u(:,1,1) = ones(ng,1);  
        
        fh_q = zeros(ng,nch,nch,2);        
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;
        fh_uh(:,4,:) = gam*gam1*M2*[-uh(:,4)./uh(:,1).^2, zero, zero, 1./uh(:,1)];
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 4 % Prescribed pressure
        %p  = uinf(:,1);        
        pinf = (uinf(:,4)-0.5)*(gam-1);
        
        uinf = [ u(:,1), u(:,2), u(:,3), pinf./(gam-1) + 0.5*u(:,2).*u(:,2)./u(:,1) + 0.5*u(:,3).*u(:,3)./u(:,1)];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zeros(ng,1)];                 
        
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An-permute(mapContractK(an-An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2]) - 2*An;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));    
    case 5 % periodic
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time); 
    case 6 % Dirichlet bc   
        %max(abs(uinf(:)))
        fh = uh-uinf;
        fh_udg = zeros(ng,nch,3*nch);
        fh_uh = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));      
    case 7
        pinf = (uinf(:,4)-0.5)*(gam-1);        
        uinf = [ u(:,1), u(:,2), u(:,3), pinf];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = 0; 

        uhinf = [uh(:,1), uh(:,2), uh(:,3), (gam-1)*(uh(:,4) - 0.5*uh(:,2).*uh(:,2)./uh(:,1) - 0.5*uh(:,3).*uh(:,3)./uh(:,1))];
        uhinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uhinfu(:,4,:) = (gam-1)*[0.5*(uh(:,2).*uh(:,2)+uh(:,3).*uh(:,3))./(uh(:,1).*uh(:,1)), -uh(:,2)./uh(:,1), -uh(:,3)./uh(:,1), ones(ng,1)];                 
        
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(An,uinf-uhinf,2,3,1,2,[],1),[2 1]);
        fh_u = permute(mapContractK(An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
        fh_uh = permute(mapContractK(Anm,uinf-uhinf,[2 4],3,1,2,[],1),[3 1 2])-permute(mapContractK(An,uhinfu,2,3,1,2,3,1),[3 1 2]);                
    case 8  % symmetry          
        un = u(:,2).*nl(:,1) + u(:,3).*nl(:,2);        
        uinf = u;
        uinf(:,2) = uinf(:,2) - nl(:,1).*un;
        uinf(:,3) = uinf(:,3) - nl(:,2).*un;
        
        [fh4,fh4_udg,fh4_uh] = fhat(nl,p,udg,uh,param,time);
        
        r    = udg(:,1);
        ru   = udg(:,2);
        rv   = udg(:,3);
        rx   = udg(:,5);
        rux  = udg(:,6);
        rvx  = udg(:,7);
        rEx  = udg(:,8);
        ry   = udg(:,9);
        ruy  = udg(:,10);
        rvy  = udg(:,11);
        rEy  = udg(:,12);
        nx   = nl(:,1);
        ny   = nl(:,2);
        r1   = 1./r;
        u    = ru.*r1;
        v    = rv.*r1;
        q    = 0.5*(u.*u+v.*v);        
        ux  = (rux - rx.*u).*r1;
        vx  = (rvx - rx.*v).*r1;
        qx  = u.*ux + v.*vx;
        px  = gam1*(rEx - rx.*q - r.*qx);        
        uy  = (ruy - ry.*u).*r1;
        vy  = (rvy - ry.*v).*r1;        
        qy  = u.*uy + v.*vy;
        py  = gam1*(rEy - ry.*q - r.*qy);             
        fh1 = px.*nx+py.*ny+(udg(:,1)-uh(:,1));
        fh1_udg = zeros(ng,1,nc);
        fh1_udg(:,1,:) = [(r.^3 - gam1.*nx.*ru.^2.*rx - gam1.*nx.*rv.^2.*rx - gam1.*ny.*ru.^2.*ry - gam1.*ny.*rv.^2.*ry + gam1.*nx.*r.*ru.*rux + gam1.*ny.*r.*ru.*ruy + gam1.*nx.*r.*rv.*rvx + gam1.*ny.*r.*rv.*rvy)./r.^3 ...
                            -(gam1.*(nx.*r.*rux + ny.*r.*ruy - nx.*ru.*rx - ny.*ru.*ry))./r.^2 ...
                            -(gam1.*(nx.*r.*rvx + ny.*r.*rvy - nx.*rv.*rx - ny.*rv.*ry))./r.^2 ...
                            0*nx ...
                            (gam1.*nx.*(ru.^2 + rv.^2))./(2.*r.^2) ...
                            -(gam1.*nx.*ru)./r ...
                            -(gam1.*nx.*rv)./r ...
                            gam1.*nx ...
                            (gam1.*ny.*(ru.^2 + rv.^2))./(2.*r.^2) ...
                            -(gam1.*ny.*ru)./r ...
                            -(gam1.*ny.*rv)./r ...
                            gam1.*ny];
        fh1_uh = zeros(ng,1,nch);
        fh1_uh(:,:,1) = -1;
                
        fh = uinf - uh;
        fh(:,1) = fh1;
        fh(:,4) = fh4(:,4);
        
        fh_u = zeros(ng,nch,nch);
        fh_u(:,2,2:3) = [ones(ng,1)-nl(:,1).*nl(:,1), -nl(:,1).*nl(:,2)]; 
        fh_u(:,3,2:3) = [ -nl(:,1).*nl(:,2), ones(ng,1)-nl(:,2).*nl(:,2)];
        fh_u(:,1,:) = fh1_udg(:,1,1:nch);
        fh_u(:,4,:) = fh4_udg(:,4,1:nch);
   
        fh_q = zeros(ng,nch,nch,2);
        fh_q(:,1,:) = fh1_udg(:,1,nch+1:3*nch);
        fh_q(:,4,:) = fh4_udg(:,4,nch+1:3*nch);
        
        fh_uh = zeros(ng,nch,nch);        
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,1,:) = fh1_uh(:,1,:);
        fh_uh(:,4,:) = fh4_uh(:,4,:);
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));                
    case 10 % Prescribing everything except pressure (taken from FM/rans)
        uinf(:,4) = u(:,4);
        uinfu = zeros(ng,nch,nch);      
        uinfu(:,4,4) = 1;
        
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = permute(mapContractK(An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
        fh_uh = permute(mapContractK(Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-An;      
    case 11
        mu = 1/Re; 
        Rex = Re;        
        Cf = 0.027/(Rex^(1/7));
        tauw = Cf/2;
        utau = sqrt(tauw);
        y = p(:,2);
        yplus = y*utau/mu;
        uplus = 0*yplus;
        ub = uinf;
        for i = 1:length(y)            
            uplus(i) = finduplus(yplus(i),1);
            ub(i,2) = uplus(i)*utau;    
            if ub(i,2)>1
                ub(i,2)=1;
            end
        end        
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);        
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,ub-uh,2,3,1,2,[],1),[2 1]);

        fh_u = an+An;
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,ub-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 12
        mu = 1/Re; 
        Rex = 2*Re;        
        Cf = 0.027/(Rex^(1/7));
        tauw = Cf/2;
        utau = sqrt(tauw);
        y = p(:,2);
        yplus = y*utau/mu;
        uplus = 0*yplus;
        ub = uinf;
        for i = 1:length(y)            
            uplus(i) = finduplus(yplus(i),1);
            ub(i,2) = uplus(i)*utau;    
            if ub(i,2)>1
                ub(i,2)=1;
            end
        end

        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);        
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,ub-uh,2,3,1,2,[],1),[2 1]);

        fh_u = an+An;
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,ub-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));     
    case 15 % Isothermal wall
        
        %one = ones(ng,1);        
        %zero = zeros(ng,1);
        
        u1 = uh(:,2)./uh(:,1);
        u2 = uh(:,3)./uh(:,1);        
        fh(:,1) =  0.978-uh(:,1);
        x = p(:,1);
        y = p(:,2);
        omega = -1+0.0*y; % = -1 + 0.5*sin(theta)        
        fh(:,2) =-omega.*y-u1;
        fh(:,3) = omega.*x-u2;        
%         fh(:,2) =  p(:,2)-u1;
%         fh(:,3) = -p(:,1)-u2;
                        
        fh_u = zeros(ng,nch,nch);
        %fh_u(:,1,1) = ones(ng,1);          
        fh_q = zeros(ng,nch,nch,2);        
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;        
        
        fh_uh(:,2,1) = uh(:,2)./(uh(:,1).^2);
        fh_uh(:,2,2) = -1./(uh(:,1));
        fh_uh(:,3,1) = uh(:,3)./(uh(:,1).^2);
        fh_uh(:,3,3) = -1./(uh(:,1));

%         pinf = (uinf(:,4)-0.5)*(gam-1);
%         fh(:,4) = (pinf-0.5)/(gam-1) - uh(:,4)./uh(:,1);
%         fh_uh(:,4,1) = uh(:,4)./uh(:,1).^2;
%         fh_uh(:,4,4) = -1./uh(:,1);
        
        pres = (gam-1)*(uh(:,4) - 0.5*(uh(:,2).*u1 + uh(:,3).*u2));
        pinf = (uinf(:,4)-0.5)*(gam-1);
        fh(:,4) = pinf - pres./uh(:,1);

        r = uh(:,1);
        ru = uh(:,2);
        rv = uh(:,3);
        rE = uh(:,4);
        fh_uh(:,4,1) = -((gam - 1)*(ru.^2./(2*r) + rv.^2./(2*r) - rE))./r.^2 - ((ru.^2./(2*r.^2) + rv.^2./(2*r.^2))*(gam - 1))./r;        
        %fh_uh(:,4,1) = -(gam-1)*0.5*(uh(:,2).^2+uh(:,3).^2)./uh(:,1).^2; %gam*gam1*M2*[-uh(:,4)./uh(:,1).^2, zero, zero, 1./uh(:,1)];
        fh_uh(:,4,2) = (gam-1)*u1./r;
        fh_uh(:,4,3) = (gam-1)*u2./r;
        fh_uh(:,4,4) = (1-gam)./r;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));        
    case 16 % Isothermal wall
        u1 = uh(:,2)./uh(:,1);
        u2 = uh(:,3)./uh(:,1);                
        x = p(:,1);
        y = p(:,2);
        r = sqrt(x.^2 + y.^2);
        vtheta = -(4-r.^2)./(3*r);
        vx = -vtheta.*y./r;
        vy = vtheta.*x./r;        
        fh(:,2) = vx-u1;
        fh(:,3) = vy-u2;        
        
        fh(:,1) = u(:,1) - uh(:,1);                
        fh(:,4) = (q(:,4,1).*nl(:,1) + q(:,4,2).*nl(:,2)); 
        fh_u = zeros(ng,nch,nch);
        fh_u(:,1,1) = ones(ng,1);                  
        fh_q = zeros(ng,nch,nch,2);        
        fh_q(:,4,4,1) = nl(:,1);
        fh_q(:,4,4,2) = nl(:,2);
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;                
        fh_uh(:,2,1) = uh(:,2)./(uh(:,1).^2);
        fh_uh(:,2,2) = -1./(uh(:,1));
        fh_uh(:,3,1) = uh(:,3)./(uh(:,1).^2);
        fh_uh(:,3,3) = -1./(uh(:,1));                  
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));                
    otherwise
        error('unknown boundary type');
end



function up = finduplus(yp,u0)

kappa = 0.41; B = 5.0; 
C = exp(-kappa*B); D = 0;
uk = kappa*u0;
y0 = u0 + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;

tol = 1e-6;
if abs(yp-y0)<tol
    up = u0;
else
    % find lower (u0) and upper (u1) bounds
    if y0<yp
        u1 = max(1,2*u0);
        uk = kappa*u1;
        y1 = u1 + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;         
        while y1<yp
            u1 = 2*u1;
            uk = kappa*u1;
            y1 = u1 + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;    
        end
    else
        u1 = u0;
        u0 = u1/2;
        uk = kappa*u0;
        y0 = u0 + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;
        while y0>yp
            u0 = u0/2;
            uk = kappa*u0;
            y0 = u0 + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;
        end
    end
    
    % bisection
    ua = 0.5*(u0+u1);        
    uk = kappa*ua;
    ya = ua + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;    
    while abs(yp-ya)>tol
        if ya<yp
            u0 = ua;
        else
            u1 = ua;
        end
        ua = 0.5*(u0+u1);
        uk = kappa*ua;
        ya = ua + C*(exp(uk) - 1 - uk - 0.5*uk.^2 - uk.^3/6) + D;    
    end    
    up = ua;
end












