function [f,f_udg] = flux2(p,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
%
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

[ng,~] = size(udg);
nch = 4;

zero = zeros(ng,1);
%one  = ones(ng,1);

gam  = param{1};
gam1 = gam - 1.0;
                                             
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);
rx   = udg(:,5);
rux  = udg(:,6);
rvx  = udg(:,7);
rEx  = udg(:,8);
ry   = udg(:,9);
ruy  = udg(:,10);
rvy  = udg(:,11);
rEy  = udg(:,12);
av = p(:,3);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
q   = 0.5*(uv.*uv+vv.*vv);
p    = gam1*(rE-r.*q);
h    = E+p.*r1;
                                        
f = zeros(ng,nch,2);
f(:,:,1) = [ru+av.*rx, ru.*uv+p+av.*rux, rv.*uv+av.*rvx,   ru.*h+av.*(rEx)];
f(:,:,2) = [rv+av.*ry, ru.*vv+av.*ruy,   rv.*vv+p+av.*rvy, rv.*h+av.*(rEy)];
% f(:,:,1) = [ru, ru.*uv+p, rv.*uv,   ru.*h];
% f(:,:,2) = [rv, ru.*vv,   rv.*vv+p, rv.*h];

if length(param) > 4
  Minf = param{2};
  Re   = param{4};
  Re1  = 1/Re;
  Pr   = param{5};
  Tref = param{6};  
  Tinf = 1/(gam*gam1*Minf^2);
  c23  = 2/3;  
    
% fc*T = gam/(gam1*Re*Pr)*p/r;  
  Ts = 110.4;
  T = p./(gam1*r);  
  Tphy =  (Tref/Tinf) * T;
  muRef = Re1;
  mu = muRef*(Tphy./Tref).^(3/2) .* (Tref + Ts)./(Tphy + Ts);  
  fc = mu*gam/Pr;
  %mu = getViscosity(muRef,Tref,Tphy,1);  
  
%   mu_r = (2^(1/2)*muRef*(5*Tref + 552)*(-(ru^2 + rv^2 - 2*r*rE)/(Tinf*r^2))^(1/2)*(ru^2 + rv^2 - r*rE)*(3312*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2))/(2*r*(1104*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2)^2); 
%   mu_ru = -(2^(1/2)*muRef*ru*(5*Tref + 552)*(-(ru^2 + rv^2 - 2*r*rE)/(Tinf*r^2))^(1/2)*(3312*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2))/(2*(1104*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2)^2);
%   mu_rv = -(2^(1/2)*muRef*rv*(5*Tref + 552)*(-(ru^2 + rv^2 - 2*r*rE)/(Tinf*r^2))^(1/2)*(3312*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2))/(2*(1104*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2)^2); 
%   mu_rE = (2^(1/2)*muRef*r*(5*Tref + 552)*(-(ru^2 + rv^2 - 2*r*rE)/(Tinf*r^2))^(1/2)*(3312*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2))/(2*(1104*Tinf*r^2 + 10*Tref*rE*r - 5*Tref*ru^2 - 5*Tref*rv^2)^2);
  
  ux  = (rux - rx.*uv).*r1;
  vx  = (rvx - rx.*vv).*r1;
  %Ex  = (rEx - rx.*E).*r1;
  qx  = uv.*ux + vv.*vx;
  px  = gam1*(rEx - rx.*q - r.*qx);
  Tx  = (1/gam1)*(px.*r - p.*rx).*r1.^2;
     
  uy  = (ruy - ry.*uv).*r1;
  vy  = (rvy - ry.*vv).*r1;
  %Ey  = (rEy - ry.*E).*r1;
  qy  = uv.*uy + vv.*vy;
  py  = gam1*(rEy - ry.*q - r.*qy);
  Ty  = (1/gam1)*(py.*r - p.*ry).*r1.^2;

  txx = c23*mu.*(2*ux - vy);
  txy = mu.*(uy + vx);
  tyy = c23*mu.*(2*vy - ux);

  fv = zeros(ng,nch,2);
  fv(:,:,1) = [zero, txx, txy, uv.*txx + vv.*txy + fc.*Tx];
  fv(:,:,2) = [zero, txy, tyy, uv.*txy + vv.*tyy + fc.*Ty];
  f = f + fv;
end

if nargout>1
    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[zeros(ng,1), 0.5*((3-gam)*uv.*uv-gam1*vv.*vv), uv.*vv, gam*E.*uv-2*gam1*uv.*q];
    f_u(:,:,1,2) = -[-ones(ng,1), (gam-3)*uv, -vv, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)];
    f_u(:,:,1,3) = -[zeros(ng,1), gam1*vv, -uv, gam1*uv.*vv];
    f_u(:,:,1,4) = -[zeros(ng,1), -gam1*ones(ng,1), zeros(ng,1), -gam*uv];    
    
    f_u(:,:,2,1) = -[zeros(ng,1), uv.*vv, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv), gam*E.*vv-2*gam1*vv.*q];
    f_u(:,:,2,2) = -[zeros(ng,1), -vv, gam1*uv, gam1*uv.*vv];
    f_u(:,:,2,3) = -[-ones(ng,1), -uv, (gam-3)*vv,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv) ];
    f_u(:,:,2,4) = -[zeros(ng,1), zeros(ng,1), -gam1*ones(ng,1), -gam*vv];

    %f_udg = f_u;    
    f_q = zeros(ng,nch,2,2*nch);
    for i = 1:nch
        f_q(:,i,1,i) = av;
        f_q(:,i,2,nch+i) = av;
    end    
        
    if length(param) > 4
    u_r  = -uv.*r1;
    u_ru =  r1;
    v_r  = -vv.*r1;
    v_rv =  r1;
    
    mu_r = (2.^(1./2).*muRef.*(5.*Tref + 552).*(-(ru.^2 + rv.^2 - 2.*r.*rE)./(Tinf.*r.^2)).^(1./2).*(ru.^2 + rv.^2 - r.*rE).*(3312.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2))./(2.*r.*(1104.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2).^2); 
    mu_ru = -(2.^(1./2).*muRef.*ru.*(5.*Tref + 552).*(-(ru.^2 + rv.^2 - 2.*r.*rE)./(Tinf.*r.^2)).^(1./2).*(3312.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2))./(2.*(1104.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2).^2);
    mu_rv = -(2.^(1./2).*muRef.*rv.*(5.*Tref + 552).*(-(ru.^2 + rv.^2 - 2.*r.*rE)./(Tinf.*r.^2)).^(1./2).*(3312.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2))./(2.*(1104.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2).^2); 
    mu_rE = (2.^(1./2).*muRef.*r.*(5.*Tref + 552).*(-(ru.^2 + rv.^2 - 2.*r.*rE)./(Tinf.*r.^2)).^(1./2).*(3312.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2))./(2.*(1104.*Tinf.*r.^2 + 10.*Tref.*rE.*r - 5.*Tref.*ru.^2 - 5.*Tref.*rv.^2).^2);
    fc_r = mu_r*gam/Pr;
    fc_ru = mu_ru*gam/Pr;
    fc_rv = mu_rv*gam/Pr;
    fc_rE = mu_rE*gam/Pr;
    
    tm = c23.*(2*ux - vy);
    txx_r  =  c23*mu.*((4*ru.*rx-2*rv.*ry)-r.*(2*rux-rvy)).*r1.^3 + mu_r.*tm;
    txx_ru = -mu.*(c23*2).*rx.*r1.^2 + mu_ru.*tm;
    txx_rv =  mu.*c23.*ry.*r1.^2 + + mu_rv.*tm;
    txx_rE =  zero + mu_rE.*tm;

    txx_rx  = -mu.*c23*2.*ru.*r1.^2;
    txx_rux =  mu.*c23*2.*r1;
    txx_rvx =  zero;
    txx_rEx =  zero;

    txx_ry  =  mu.*c23.*rv.*r1.^2;
    txx_ruy =  zero;
    txx_rvy = -mu.*c23.*r1;
    txx_rEy =  zero;

    tm = (uy + vx);
    txy_r  =  mu.*(2*(ru.*ry+rv.*rx)-r.*(ruy+rvx)).*r1.^3 + mu_r.*tm;
    txy_ru = -mu.*ry.*r1.^2 + mu_ru.*tm;
    txy_rv = -mu.*rx.*r1.^2 + mu_rv.*tm;
    txy_rE =  zero + mu_rE.*tm;

    txy_rx  = -mu.*rv.*r1.^2;
    txy_rux =  zero;
    txy_rvx =  mu.*r1;
    txy_rEx =  zero;

    txy_ry  = -mu.*ru.*r1.^2;
    txy_ruy =  mu.*r1;
    txy_rvy =  zero;
    txy_rEy =  zero;

    tm = c23.*(2*vy - ux);
    tyy_r  =  mu.*c23.*((4*rv.*ry-2*ru.*rx)-r.*(2*rvy-rux)).*r1.^3 + mu_r.*tm;
    tyy_ru =  mu.*c23.*rx.*r1.^2 + mu_ru.*tm;
    tyy_rv = -mu.*(c23*2).*ry.*r1.^2 + mu_rv.*tm;
    tyy_rE =  zero + mu_rE.*tm;

    tyy_rx  =  mu.*c23.*ru.*r1.^2;
    tyy_rux = -mu.*c23.*r1;
    tyy_rvx =  zero;
    tyy_rEx =  zero;

    tyy_ry  = -mu.*(c23*2).*rv.*r1.^2;
    tyy_ruy =  zero;
    tyy_rvy =  mu.*(c23*2).*r1;
    tyy_rEy =  zero;

    Tx_r  = -(rEx.*r.^2-2*rux.*r.*ru-2*rvx.*r.*rv-2*rE.*rx.*r+3*rx.*(ru.^2+rv.^2)).*r1.^4;
    Tx_ru = -(r.*rux-2*ru.*rx).*r1.^3;
    Tx_rv = -(r.*rvx-2*rv.*rx).*r1.^3;
    Tx_rE = -rx.*r1.^2;

    Tx_rx  =  (ru.^2+rv.^2-r.*rE).*r1.^3;
    Tx_rux = -ru.*r1.^2;
    Tx_rvx = -rv.*r1.^2;
    Tx_rEx =  r1;

    Ty_r  = -(rEy.*r.^2-2*ruy.*r.*ru-2*rvy.*r.*rv-2*rE.*ry.*r+3*ry.*(ru.^2+rv.^2)).*r1.^4;
    Ty_ru = -(r.*ruy-2*ru.*ry).*r1.^3;
    Ty_rv = -(r.*rvy-2*rv.*ry).*r1.^3;
    Ty_rE = -ry.*r1.^2;

    Ty_ry  = Tx_rx;
    Ty_ruy = Tx_rux;
    Ty_rvy = Tx_rvx;
    Ty_rEy = Tx_rEx;

    fv_u = zeros(ng,nch,2,nch);
    fv_u(:,:,1,1) = [zero, txx_r , txy_r , u_r.*txx  + uv.*txx_r  + v_r.*txy  + vv.*txy_r  + fc.*Tx_r + fc_r.*Tx];
    fv_u(:,:,1,2) = [zero, txx_ru, txy_ru, u_ru.*txx + uv.*txx_ru             + vv.*txy_ru + fc.*Tx_ru + fc_ru.*Tx];
    fv_u(:,:,1,3) = [zero, txx_rv, txy_rv,             uv.*txx_rv + v_rv.*txy + vv.*txy_rv + fc.*Tx_rv + fc_rv.*Tx];
    fv_u(:,:,1,4) = [zero, txx_rE, txy_rE,             uv.*txx_rE             + vv.*txy_rE + fc.*Tx_rE + fc_rE.*Tx];

    fv_u(:,:,2,1) = [zero, txy_r , tyy_r , u_r.*txy  + uv.*txy_r  + v_r.*tyy  + vv.*tyy_r  + fc.*Ty_r + fc_r.*Ty];
    fv_u(:,:,2,2) = [zero, txy_ru, tyy_ru, u_ru.*txy + uv.*txy_ru             + vv.*tyy_ru + fc.*Ty_ru + fc_ru.*Ty];
    fv_u(:,:,2,3) = [zero, txy_rv, tyy_rv,             uv.*txy_rv + v_rv.*tyy + vv.*tyy_rv + fc.*Ty_rv + fc_rv.*Ty];
    fv_u(:,:,2,4) = [zero, txy_rE, tyy_rE,             uv.*txy_rE             + vv.*tyy_rE + fc.*Ty_rE + fc_rE.*Ty];

    fv_q = zeros(ng,nch,2,nch,2);
    fv_q(:,:,1,1,1) = [zero, txx_rx , txy_rx , uv.*txx_rx  + vv.*txy_rx  + fc.*Tx_rx ];
    fv_q(:,:,1,2,1) = [zero, txx_rux, txy_rux, uv.*txx_rux + vv.*txy_rux + fc.*Tx_rux];
    fv_q(:,:,1,3,1) = [zero, txx_rvx, txy_rvx, uv.*txx_rvx + vv.*txy_rvx + fc.*Tx_rvx];
    fv_q(:,:,1,4,1) = [zero, txx_rEx, txy_rEx, uv.*txx_rEx + vv.*txy_rEx + fc.*Tx_rEx];
    fv_q(:,:,1,1,2) = [zero, txx_ry , txy_ry , uv.*txx_ry  + vv.*txy_ry             ];
    fv_q(:,:,1,2,2) = [zero, txx_ruy, txy_ruy, uv.*txx_ruy + vv.*txy_ruy            ];
    fv_q(:,:,1,3,2) = [zero, txx_rvy, txy_rvy, uv.*txx_rvy + vv.*txy_rvy            ];
    fv_q(:,:,1,4,2) = [zero, txx_rEy, txy_rEy, uv.*txx_rEy + vv.*txy_rEy            ];

    fv_q(:,:,2,1,1) = [zero, txy_rx , tyy_rx , uv.*txy_rx  + vv.*tyy_rx             ];
    fv_q(:,:,2,2,1) = [zero, txy_rux, tyy_rux, uv.*txy_rux + vv.*tyy_rux            ];
    fv_q(:,:,2,3,1) = [zero, txy_rvx, tyy_rvx, uv.*txy_rvx + vv.*tyy_rvx            ];
    fv_q(:,:,2,4,1) = [zero, txy_rEx, tyy_rEx, uv.*txy_rEx + vv.*tyy_rEx            ];
    fv_q(:,:,2,1,2) = [zero, txy_ry , tyy_ry , uv.*txy_ry  + vv.*tyy_ry  + fc.*Ty_ry ];
    fv_q(:,:,2,2,2) = [zero, txy_ruy, tyy_ruy, uv.*txy_ruy + vv.*tyy_ruy + fc.*Ty_ruy];
    fv_q(:,:,2,3,2) = [zero, txy_rvy, tyy_rvy, uv.*txy_rvy + vv.*tyy_rvy + fc.*Ty_rvy];
    fv_q(:,:,2,4,2) = [zero, txy_rEy, tyy_rEy, uv.*txy_rEy + vv.*tyy_rEy + fc.*Ty_rEy];
    
    fv_q = reshape(fv_q, size(f_q));
    f_u = f_u + fv_u;
    f_q = f_q + fv_q;
    end
    
    f_udg = cat(4,f_u,f_q);
end
