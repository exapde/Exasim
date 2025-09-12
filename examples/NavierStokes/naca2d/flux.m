function [f,f_udg] = flux(p,udg,param,time)
% FLUX Volume flux function
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

nc = size(udg,2);

if nc == 12
    [f,f_udg] = flux2d(p,udg,param,time);
elseif nc == 20
    [f,f_udg] = flux3d(p,udg,param,time);
end

% [ng,nc] = size(udg);
% nch = nc/3;
% zero = zeros(ng,1);
% one  = ones(ng,1);
% 
% gam  = param{1};
% gam1 = gam-1.0;
% Re   = param{3};
% Re1  = 1/Re;
% Pr   = param{4};
% Minf = param{5};
% M2   = Minf^2;
% 
% c23  = 2/3;
% fc = 1/(gam1*M2*Re*Pr);
%                                              
% r    = udg(:,1);
% ru   = udg(:,2);
% rv   = udg(:,3);
% rE   = udg(:,4);
% 
% rx   = udg(:,5);
% rux  = udg(:,6);
% rvx  = udg(:,7);
% rEx  = udg(:,8);
% 
% ry   = udg(:,9);
% ruy  = udg(:,10);
% rvy  = udg(:,11);
% rEy  = udg(:,12);
% 
% r1   = 1./r;
% u    = ru.*r1;
% v    = rv.*r1;
% E    = rE.*r1;
% q    = 0.5*(u.*u+v.*v);
% p    = gam1*(rE-r.*q);
% h    = E+p.*r1;
%                                         
% f = zeros(ng,nch,2);
% f(:,:,1) = [ru, ru.*u+p, rv.*u,   ru.*h];
% f(:,:,2) = [rv, ru.*v,   rv.*v+p, rv.*h];
% 
% f_u = zeros(ng,nch,2,nch);
% f_u(:,:,1,1) = [zero, 0.5*((gam-3)*u.*u+gam1*v.*v), -u.*v,         2*gam1*u.*q-gam*E.*u];
% f_u(:,:,1,2) = [ one,                    (3-gam)*u,     v, gam*E-0.5*gam1*(3*u.*u+v.*v)];
% f_u(:,:,1,3) = [zero,                      -gam1*v,     u,                   -gam1*u.*v];
% f_u(:,:,1,4) = [zero,                     gam1*one,  zero,                        gam*u];
% 
% f_u(:,:,2,1) = [zero, -u.*v, 0.5*((gam-3)*v.*v+gam1*u.*u),         2*gam1*v.*q-gam*E.*v];
% f_u(:,:,2,2) = [zero,     v,                      -gam1*u,                   -gam1*u.*v];
% f_u(:,:,2,3) = [ one,     u,                    (3-gam)*v, gam*E-0.5*gam1*(3*v.*v+u.*u)];
% f_u(:,:,2,4) = [zero,  zero,                     gam1*one,                        gam*v];
% 
% u_r  = -u.*r1;
% u_ru =  r1;
% v_r  = -v.*r1;
% v_rv =  r1;
% 
% ux  = (rux - rx.*u).*r1;
% vx  = (rvx - rx.*v).*r1;
% Ex  = (rEx - rx.*E).*r1;
% qx  = u.*ux + v.*vx;
% px  = gam1*(rEx - rx.*q - r.*qx);
% Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;
% 
% uy  = (ruy - ry.*u).*r1;
% vy  = (rvy - ry.*v).*r1;
% Ey  = (rEy - ry.*E).*r1;
% qy  = u.*uy + v.*vy;
% py  = gam1*(rEy - ry.*q - r.*qy);
% Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;
% 
% txx = Re1*c23*(2*ux - vy);
% txy = Re1*(uy + vx);
% tyy = Re1*c23*(2*vy - ux);
% 
% % Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;
% % Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;
% 
% fv = zeros(ng,nch,2);
% fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
% fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];
% 
% txx_r  =  Re1*c23*((4*ru.*rx-2*rv.*ry)-r.*(2*rux-rvy)).*r1.^3;
% txx_ru = -Re1*c23*2*rx.*r1.^2;
% txx_rv =  Re1*c23*ry.*r1.^2;
% txx_rE =  zero;
% 
% txx_rx  = -Re1*c23*2*ru.*r1.^2;
% txx_rux =  Re1*c23*2*r1;
% txx_rvx =  zero;
% txx_rEx =  zero;
% 
% txx_ry  =  Re1*c23*rv.*r1.^2;
% txx_ruy =  zero;
% txx_rvy = -Re1*c23.*r1;
% txx_rEy =  zero;
% 
% txy_r  =  Re1*(2*(ru.*ry+rv.*rx)-r.*(ruy+rvx)).*r1.^3;
% txy_ru = -Re1*ry.*r1.^2;
% txy_rv = -Re1*rx.*r1.^2;
% txy_rE =  zero;
% 
% txy_rx  = -Re1*rv.*r1.^2;
% txy_rux =  zero;
% txy_rvx =  Re1*r1;
% txy_rEx =  zero;
% 
% txy_ry  = -Re1*ru.*r1.^2;
% txy_ruy =  Re1*r1;
% txy_rvy =  zero;
% txy_rEy =  zero;
% 
% tyy_r  =  Re1*c23*((4*rv.*ry-2*ru.*rx)-r.*(2*rvy-rux)).*r1.^3;
% tyy_ru =  Re1*c23*rx.*r1.^2;
% tyy_rv = -Re1*c23*2*ry.*r1.^2;
% tyy_rE =  zero;
% 
% tyy_rx  =  Re1*c23*ru.*r1.^2;
% tyy_rux = -Re1*c23*r1;
% tyy_rvx =  zero;
% tyy_rEx =  zero;
% 
% tyy_ry  = -Re1*c23*2*rv.*r1.^2;
% tyy_ruy =  zero;
% tyy_rvy =  Re1*c23*2*r1;
% tyy_rEy =  zero;
% 
% Tx_r  = -M2*gam*gam1*(rEx.*r.^2-2*rux.*r.*ru-2*rvx.*r.*rv-2*rE.*rx.*r+3*rx.*(ru.^2+rv.^2)).*r1.^4;
% Tx_ru = -M2*gam*gam1*(r.*rux-2*ru.*rx).*r1.^3;
% Tx_rv = -M2*gam*gam1*(r.*rvx-2*rv.*rx).*r1.^3;
% Tx_rE = -M2*gam*gam1*rx.*r1.^2;
% 
% Tx_rx  =  M2*gam*gam1*(ru.^2+rv.^2-r.*rE).*r1.^3;
% Tx_rux = -M2*gam*gam1*ru.*r1.^2;
% Tx_rvx = -M2*gam*gam1*rv.*r1.^2;
% Tx_rEx =  M2*gam*gam1*r1;
% 
% Ty_r  = -M2*gam*gam1*(rEy.*r.^2-2*ruy.*r.*ru-2*rvy.*r.*rv-2*rE.*ry.*r+3*ry.*(ru.^2+rv.^2)).*r1.^4;
% Ty_ru = -M2*gam*gam1*(r.*ruy-2*ru.*ry).*r1.^3;
% Ty_rv = -M2*gam*gam1*(r.*rvy-2*rv.*ry).*r1.^3;
% Ty_rE = -M2*gam*gam1*ry.*r1.^2;
% 
% Ty_ry  = Tx_rx;
% Ty_ruy = Tx_rux;
% Ty_rvy = Tx_rvx;
% Ty_rEy = Tx_rEx;
% 
% % fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
% % fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];
% 
% fv_u = zeros(ng,nch,2,nch);
% fv_u(:,:,1,1) = [zero, txx_r , txy_r , u_r.*txx  + u.*txx_r  + v_r.*txy  + v.*txy_r  + fc*Tx_r ];
% fv_u(:,:,1,2) = [zero, txx_ru, txy_ru, u_ru.*txx + u.*txx_ru             + v.*txy_ru + fc*Tx_ru];
% fv_u(:,:,1,3) = [zero, txx_rv, txy_rv,             u.*txx_rv + v_rv.*txy + v.*txy_rv + fc*Tx_rv];
% fv_u(:,:,1,4) = [zero, txx_rE, txy_rE,             u.*txx_rE             + v.*txy_rE + fc*Tx_rE];
% 
% fv_u(:,:,2,1) = [zero, txy_r , tyy_r , u_r.*txy  + u.*txy_r  + v_r.*tyy  + v.*tyy_r  + fc*Ty_r];
% fv_u(:,:,2,2) = [zero, txy_ru, tyy_ru, u_ru.*txy + u.*txy_ru             + v.*tyy_ru + fc*Ty_ru];
% fv_u(:,:,2,3) = [zero, txy_rv, tyy_rv,             u.*txy_rv + v_rv.*tyy + v.*tyy_rv + fc*Ty_rv];
% fv_u(:,:,2,4) = [zero, txy_rE, tyy_rE,             u.*txy_rE             + v.*tyy_rE + fc*Ty_rE];
% 
% fv_q = zeros(ng,nch,2,nch,2);
% fv_q(:,:,1,1,1) = [zero, txx_rx , txy_rx , u.*txx_rx  + v.*txy_rx  + fc*Tx_rx ];
% fv_q(:,:,1,2,1) = [zero, txx_rux, txy_rux, u.*txx_rux + v.*txy_rux + fc*Tx_rux];
% fv_q(:,:,1,3,1) = [zero, txx_rvx, txy_rvx, u.*txx_rvx + v.*txy_rvx + fc*Tx_rvx];
% fv_q(:,:,1,4,1) = [zero, txx_rEx, txy_rEx, u.*txx_rEx + v.*txy_rEx + fc*Tx_rEx];
% fv_q(:,:,1,1,2) = [zero, txx_ry , txy_ry , u.*txx_ry  + v.*txy_ry             ];
% fv_q(:,:,1,2,2) = [zero, txx_ruy, txy_ruy, u.*txx_ruy + v.*txy_ruy            ];
% fv_q(:,:,1,3,2) = [zero, txx_rvy, txy_rvy, u.*txx_rvy + v.*txy_rvy            ];
% fv_q(:,:,1,4,2) = [zero, txx_rEy, txy_rEy, u.*txx_rEy + v.*txy_rEy            ];
% 
% fv_q(:,:,2,1,1) = [zero, txy_rx , tyy_rx , u.*txy_rx  + v.*tyy_rx             ];
% fv_q(:,:,2,2,1) = [zero, txy_rux, tyy_rux, u.*txy_rux + v.*tyy_rux            ];
% fv_q(:,:,2,3,1) = [zero, txy_rvx, tyy_rvx, u.*txy_rvx + v.*tyy_rvx            ];
% fv_q(:,:,2,4,1) = [zero, txy_rEx, tyy_rEx, u.*txy_rEx + v.*tyy_rEx            ];
% fv_q(:,:,2,1,2) = [zero, txy_ry , tyy_ry , u.*txy_ry  + v.*tyy_ry  + fc*Ty_ry ];
% fv_q(:,:,2,2,2) = [zero, txy_ruy, tyy_ruy, u.*txy_ruy + v.*tyy_ruy + fc*Ty_ruy];
% fv_q(:,:,2,3,2) = [zero, txy_rvy, tyy_rvy, u.*txy_rvy + v.*tyy_rvy + fc*Ty_rvy];
% fv_q(:,:,2,4,2) = [zero, txy_rEy, tyy_rEy, u.*txy_rEy + v.*tyy_rEy + fc*Ty_rEy];
% 
% f = f+fv;
% f_udg = cat(4,f_u+fv_u,reshape(fv_q,ng,nch,2,2*nch));