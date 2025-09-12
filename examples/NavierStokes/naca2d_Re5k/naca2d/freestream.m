nch = 4;
nd = 2;

syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms uinf1 uinf2 uinf3 uinf4 uinf5
syms uh1 uh2 uh3 uh4 uh5
syms x1 x2 x3 n1 n2 n3 av
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10
syms zero one 

param = [param1 param2 param3 param4 param5 param6 param7 param8];
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];        
    uinf = [uinf1 uinf2 uinf3 uinf4];
    uh = [uh1 uh2 uh3 uh4];
    nl = [n1 n2];
    pg = [x1 x2 x3];       
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];        
    uinf = [uinf1 uinf2 uinf3 uinf4 uinf5];
    uh = [uh1 uh2 uh3 uh4 uh5];
    nl = [n1 n2 n3];
    pg = [x1 x2 x3];      
end

gam = param(1);
gam1 = gam - 1.0;

u = udg(1:4);
r = uh(1);
ru = uh(2);
rv = uh(3);
rE = uh(4);
nx = nl(1);
ny = nl(2);

r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
h = E+p*r1;
a = sqrt(gam*p*r1);

run = ru*nx + rv*ny;
rut = -ru*ny + rv*nx;
un = run/r;
ut = rut/r;

K = [ 1 , 1 , 0 , 1 ;...
      un-a , un , 0 , un+a ;...
      ut , ut , 1 , ut ;...
      h - un*a , (1/2)*(un^2 + ut^2) , ut , h+un*a ];
Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                        -2*h + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                        -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                        h - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
T = [ 1 , 0 , 0 , 0;...
      0 , nx , ny , 0;...
      0 , -ny , nx , 0;...
      0 , 0 , 0 , 1];
Tinv = [ 1 , 0 , 0 , 0;...
         0 , nx ,-ny , 0;...
         0 , ny , nx , 0;...
         0 , 0 , 0 , 1];
Lambda = [ tanh(1e2*(un-a)) , 0 , 0 , 0 ;...
                 0 , tanh(1e2*(un)) , 0 , 0 ;...
                 0 , 0 , tanh(1e2*(un)) , 0 ;...
                 0 , 0 , 0 , tanh(1e2*(un+a)) ];
E = (K * Lambda * Kinv);
%An = (Tinv * E * T);
L = simplify(Tinv * K);
R = simplify(Kinv * T);
An = simplify(L * Lambda * R);

% An_uh   = jacobian(An(:),uh);
% E_uh   = jacobian(E(:),uh);
% matlabFunction(E(:),E_uh(:),'file','tmp.m','vars',{nl,uh,param},'outputs', {'E','E_uh'});

% freestream boundary condition
clear f fh s fb;
fb = 0.5*((u(:)+uinf(:)) + An*(u(:)-uinf(:))) - uh(:);     

matlabFunction(fb,'file','tmp.m','vars',{nl,udg,uh,uinf,param},'outputs', {'fb'});

fb_u   = jacobian(fb(:),u);
matlabFunction(fb_u(:),'file','tmp2.m','vars',{nl,udg,uh,uinf,param},'outputs', {'fb_u'});

fb_uh   = jacobian(fb(:),uh);
matlabFunction(fb_uh(:),'file','tmp3.m','vars',{nl,udg,uh,uinf,param},'outputs', {'fb_uh'});



% filename4 = "freestreambc.m";
% matlabcodegen;


