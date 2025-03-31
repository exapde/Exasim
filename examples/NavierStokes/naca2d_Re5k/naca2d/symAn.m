function [An,Anm] = symAn(nl,uh,param,flux_type)
gamma = param(1);
dim = length(nl);
if dim==2
    r = uh(1);
    ru = uh(2);
    rv = uh(3);
    rE = uh(4);
    nx = nl(1);
    ny = nl(2);
    gam1 = gamma-1;
    run = ru*nx + rv*ny;
    rut = -ru*ny + rv*nx;
    un = run/r;
    ut = rut/r;
    P = gam1*(rE - r*(1/2)*(un^2+ut^2));
    a = sqrt(gamma*P/r);
    H = rE/r + P/r;
    % A = [0 , 1 , 0 , 0 ;...
    % -un^2+(1/2)*gam1*(un^2+ut^2) , (3-gamma)*un , -gam1*ut , gam1 ;...
    % -un*ut , ut , un , 0 ;...
    % un*((1/2)*gam1*(un^2+ut^2)-H) , H-gam1*(un^2) , -gam1*un*ut , gamma*un ];
    K = [ 1 , 1 , 0 , 1 ;...
          un-a , un , 0 , un+a ;...
          ut , ut , 1 , ut ;...
          H - un*a , (1/2)*(un^2 + ut^2) , ut , H+un*a ];
    Kinv = (gam1/(2*a^2))*[ H + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                            -2*H + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                            -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                            H - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
    T = [ 1 , 0 , 0 , 0;...
          0 , nx , ny , 0;...
          0 , -ny , nx , 0;...
          0 , 0 , 0 , 1];
    Tinv = [ 1 , 0 , 0 , 0;...
             0 , nx ,-ny , 0;...
             0 , ny , nx , 0;...
             0 , 0 , 0 , 1];
    switch flux_type
        case 2 %%% Lax Friedrich
            An = [ abs(un)+a , 0 , 0 , 0 ;...
                    0 , abs(un)+a , 0 , 0 ;...
                    0 , 0 , abs(un)+a , 0 ;...
                    0 , 0 , 0 , abs(un)+a ];
            An = simplify(An);
        case 0 %%% Roe w/o entropy fix
            Lambda = [ un-a , 0 , 0 , 0 ;...
                        0 , un , 0 , 0 ;...
                        0 , 0 , un , 0 ;...
                        0 , 0 , 0 , un+a ];
            E = simplify(K * Lambda * Kinv);
            An = simplify(Tinv * E * T);
% r1 = 1/r;
% uv = ru*r1;
% vv = rv*r1;
% E = rE*r1;
% af = 0.5*(uv*uv+vv*vv);
% p = gam1*(rE-r*af);
% h = E+p*r1;
% F = [ru, ru*uv+p, rv*uv, ru*h, ...
% rv, ru*vv, rv*vv+p, rv*h];
% F = reshape(F,4,2);
% Fn = F(:,1)*nx + F(:,2)*ny;
% Bn = simplify(jacobian(Fn,uh));
% simplify(An-Bn)
        case 1 %%% Absolute_Roe w/o entropy fix
            Lambda = [ abs(un-a) , 0 , 0 , 0 ;...
                             0 , abs(un) , 0 , 0 ;...
                             0 , 0 , abs(un) , 0 ;...
                             0 , 0 , 0 , abs(un+a) ];
            E = simplify(K * Lambda * Kinv);
            An = simplify(Tinv * E * T);
        case 3
            Lambda1 = [ abs(un-a) , 0 , 0 , 0 ;...
                             0 , abs(un) , 0 , 0 ;...
                             0 , 0 , abs(un) , 0 ;...
                             0 , 0 , 0 , abs(un+a) ];
            Lambda2 = [ un-a , 0 , 0 , 0 ;...
                        0 , un , 0 , 0 ;...
                        0 , 0 , un , 0 ;...
                        0 , 0 , 0 , un+a ];
            Lambda = Lambda1+Lambda2;
            E = simplify(K * Lambda * Kinv);
            An = simplify(Tinv * E * T);
        case 4
            Lambda1 = [ abs(un-a) , 0 , 0 , 0 ;...
                             0 , abs(un) , 0 , 0 ;...
                             0 , 0 , abs(un) , 0 ;...
                             0 , 0 , 0 , abs(un+a) ];
            Lambda2 = [ un-a , 0 , 0 , 0 ;...
                        0 , un , 0 , 0 ;...
                        0 , 0 , un , 0 ;...
                        0 , 0 , 0 , un+a ];
            Lambda = Lambda1-Lambda2;
            E = simplify(K * Lambda * Kinv);
            An = simplify(Tinv * E * T);
    end
elseif dim==3
    r = uh(1);
    ru = uh(2);
    rv = uh(3);
    rw = uh(4);
    rE = uh(5);
    nx = nl(1);
    ny = nl(2);
    nz = nl(3);
    gam1 = gamma-1;
    cy = sqrt(nx^2+ny^2);
    sy = nz;
    cz = nx/cy;
    sz = ny/cy;
    Tz = [1, 0, 0, 0, 0;...
          0, cz, sz, 0, 0;...
          0, -sz, cz, 0, 0;...
          0, 0, 0, 1, 0;
          0, 0, 0, 0, 1];
    Ty = [1, 0, 0, 0, 0;...
          0, cy, 0, sy, 0;...
          0, 0, 1, 0, 0;...
          0, -sy, 0, cy, 0;
          0, 0, 0, 0, 1];
    T = Ty*Tz;
    Tinv = (Tz.')*(Ty.');
    U = [r, ru, rv, rw, rE].';
    TU = T*U;
    run = TU(2);
    rvn = TU(3);
    rwn = TU(4);
    un = run/r;
    vn = rvn/r;
    wn = rwn/r;
    P = gam1*(rE - r*(1/2)*(un^2 + vn^2 + wn^2));
    a = sqrt(gamma*P/r);
    H = rE/r + P/r;
% A = [ 0 , 1 , 0 , 0, 0 ;...
% gam1*H-un^2-a^2, (3-gamma)*un , -gam1*vn , -gam1*wn, gam1 ;...
% -un*vn , vn , un , 0 , 0 ;...
% -un*wn, wn, 0, un, 0;...
% un*((1/2)*gam1*(un^2+vn^2+wn^2)-H), H-gam1*(un^2), -gam1*un*vn , -gam1*un*wn, gamma*un];
    K = [ 1 , 1 , 0 , 0, 1 ;...
          un-a , un , 0 , 0, un+a ;...
          vn , vn , 1 , 0, vn ;...
          wn , wn , 0 , 1, wn ;...
          H - un*a , (1/2)*(un^2+vn^2+wn^2) , vn , wn, H+un*a ];
    Kinv = (gam1/(2*a^2))*[ H + (a/gam1)*(un-a) , -(un+a/gam1) , -vn , -wn, 1 ;...
                            -2*H + (4/gam1)*a^2 , 2*un , 2*vn , 2*wn, -2 ;...
                            -2*(vn*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0, 0 ;...
                            -2*wn*a^2/gam1 , 0 , 0 , 2*a^2/gam1, 0
                            H - a*(un+a)/gam1 , -un+a/gam1 , -vn , -wn, 1 ];
    if flux_type==0
        Lambda = [ un-a , 0 , 0 , 0 , 0;...
                    0 , un , 0 , 0 , 0;...
                    0 , 0 , un , 0 , 0;...
                    0 , 0 , 0 , un , 0; ...
                    0 , 0 , 0 , 0 , un+a ];
        E = simplify(K * Lambda * Kinv);
        An = simplify(Tinv * E * T);
    elseif flux_type==1
        Lambda = [ abs(un-a) , 0 , 0 , 0 , 0;...
                    0 , abs(un) , 0 , 0 , 0;...
                    0 , 0 , abs(un) , 0 , 0;...
                    0 , 0 , 0 , abs(un) , 0; ...
                    0 , 0 , 0 , 0 , abs(un+a) ];
        E = (K * Lambda * Kinv);
        An = (Tinv * E * T);
    elseif flux_type==2
        An = [ abs(un)+a , 0 , 0 , 0 , 0;...
                0 , abs(un)+a , 0 , 0 , 0;...
                0 , 0 , abs(un)+a , 0 , 0;...
                0 , 0 , 0 , abs(un)+a , 0; ...
                0 , 0 , 0 , 0 , abs(un)+a ];
    elseif flux_type==3
        Lambda1 = [ abs(un-a) , 0 , 0 , 0 , 0;...
                    0 , abs(un) , 0 , 0 , 0;...
                    0 , 0 , abs(un) , 0 , 0;...
                    0 , 0 , 0 , abs(un) , 0; ...
                    0 , 0 , 0 , 0 , abs(un+a) ];
        Lambda2 = [ un-a , 0 , 0 , 0 , 0;...
                    0 , un , 0 , 0 , 0;...
                    0 , 0 , un , 0 , 0;...
                    0 , 0 , 0 , un , 0; ...
                    0 , 0 , 0 , 0 , un+a ];
        Lambda = Lambda1+Lambda2;
        E = simplify(K * Lambda * Kinv);
        An = simplify(Tinv * E * T);
    elseif flux_type==4
        Lambda1 = [ abs(un-a) , 0 , 0 , 0 , 0;...
                    0 , abs(un) , 0 , 0 , 0;...
                    0 , 0 , abs(un) , 0 , 0;...
                    0 , 0 , 0 , abs(un) , 0; ...
                    0 , 0 , 0 , 0 , abs(un+a) ];
        Lambda2 = [ un-a , 0 , 0 , 0 , 0;...
                    0 , un , 0 , 0 , 0;...
                    0 , 0 , un , 0 , 0;...
                    0 , 0 , 0 , un , 0; ...
                    0 , 0 , 0 , 0 , un+a ];
        Lambda = Lambda1-Lambda2;
        E = simplify(K * Lambda * Kinv);
        An = simplify(Tinv * E * T);
    end
% Q = simplify(E-A)
% simplify(K*Kinv)
% simplify(Kinv*K)
% r1 = 1/r;
% uv = ru*r1;
% vv = rv*r1;
% wv = rw*r1;
% E = rE*r1;
% af = 0.5*(uv*uv+vv*vv+wv*wv);
% p = gam1*(rE-r*af);
% h = E+p*r1;
%
% F = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
% rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
% rw, ru*wv, rv*wv, rw*wv+p, rw*h];
%
% F = reshape(F,5,3);
% Fn = F(:,1)*nx + F(:,2)*ny + F(:,3)*nz;
% Bn = simplify(jacobian(Fn,uh));
% tm = simplify(An-Bn)
% An
% Bn
% a = rand/2;
% b = rand/2;
% w=subs(tm,{r,ru,rv,rw,rE,nx,ny,nz,gamma},{rand,rand,rand,rand,20-rand,a,b,sqrt(1-a^2-b^2),1.4});
% round(w)
end
syms zero one
%%% compute Jacobian
Anm = jacobian(An(:),uh);
%%% And patch with vector zero or one to have the right sizes
for ii = 1:size(An,1)
    for jj = 1:size(An,2)
        temp = ismember(symvar(An(ii,jj)),uh);
        if An(ii,jj)==0, An(ii,jj) = zero; %%% No dependency on anything
        elseif isempty(temp) || sum(temp)==0, An(ii,jj) = (An(ii,jj))*one; %%% No dependency on state vars
        end
    end
end
for ii = 1:size(Anm,1)
    for jj = 1:size(Anm,2)
        temp = ismember(symvar(Anm(ii,jj)),uh);
        if Anm(ii,jj)==0, Anm(ii,jj) = zero; %%% No dependency on anything
        elseif isempty(temp) || sum(temp)==0, Anm(ii,jj) = (Anm(ii,jj))*one; %%% No dependency on state vars
        end
    end
end
nd = length(nl);
nch = length(uh);
filename1 = ['getan' num2str(nd) 'd4' '.m'];
% generate a temporary matlab file
matlabFunction(An(:),'file','tmp.m',...
    'vars',{nl,uh,param,[zero one]},'outputs', {'An'});
% open the file and modify it
fid = fopen('tmp.m','r');
gid = fopen(filename1,'wt');
tline = fgetl(fid);
i=1;
while ischar(tline)
    str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
    str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
    str = strrep(str, 'in1', 'nl');
    str = strrep(str, 'IN1', 'NL');
    str = strrep(str, 'in2', 'uh');
    str = strrep(str, 'IN2', 'UH');
    str = strrep(str, 'in3', 'param');
    str = strrep(str, 'IN3', 'PARAM');
    str = strrep(str, ',in5)', ')');
    str = strrep(str, ',IN5)', ')');
    str = strrep(str, 'param(:,1)', 'param{1}');
    str = strrep(str, 'param(:,2)', 'param{2}');
    str = strrep(str, 'param(:,3)', 'param{3}');
    str = strrep(str, 'param(:,4)', 'param{4}');
    str = strrep(str, 'param(:,5)', 'param{5}');
    str = strrep(str, 'param(:,6)', 'param{6}');
    str = strrep(str, 'param(:,7)', 'param{7}');
    str = strrep(str, 'param(:,8)', 'param{8}');
    str = strrep(str, 'param(:,9)', 'param{9}');
    str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
    str = strrep(str, 'in5(:,2)', 'ones(ng,1)');
    if i==7
        str = '[ng,nc] = size(uh);';
        fprintf(gid, '%s\n', str);
        str = ['nch = ' num2str(nch) ';'];
        fprintf(gid, '%s\n', str);
        str = ['nd = ' num2str(nd) ';'];
    end
    fprintf(gid, '%s\n', str);
    tline = fgetl(fid);
    i=i+1;
    %disp(str)
end
str = 'An = reshape(An,ng,nch,nch);';
fprintf(gid, '%s\n', str);
str = 'Anm = reshape(Anm,ng,nch,nch,nch);';
fprintf(gid, '%s\n', str);
fclose(fid);
fclose(gid);
delete('tmp.m');

% filename1 = ['getan' num2str(nd) 'd4' '.m'];
% % generate a temporary matlab file
% matlabFunction(An(:),Anm(:),'file','tmp.m',...
%     'vars',{nl,uh,param,[zero one]},'outputs', {'f','f_udg'});
% % open the file and modify it
% fid = fopen('tmp.m','r');
% gid = fopen(filename1,'wt');
% tline = fgetl(fid);
% i=1;
% while ischar(tline)
%     str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
%     str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
%     str = strrep(str, 'in1', 'nl');
%     str = strrep(str, 'IN1', 'NL');
%     str = strrep(str, 'in2', 'uh');
%     str = strrep(str, 'IN2', 'UH');
%     str = strrep(str, 'in3', 'param');
%     str = strrep(str, 'IN3', 'PARAM');
%     str = strrep(str, ',in5)', ')');
%     str = strrep(str, ',IN5)', ')');
%     str = strrep(str, 'param(:,1)', 'param{1}');
%     str = strrep(str, 'param(:,2)', 'param{2}');
%     str = strrep(str, 'param(:,3)', 'param{3}');
%     str = strrep(str, 'param(:,4)', 'param{4}');
%     str = strrep(str, 'param(:,5)', 'param{5}');
%     str = strrep(str, 'param(:,6)', 'param{6}');
%     str = strrep(str, 'param(:,7)', 'param{7}');
%     str = strrep(str, 'param(:,8)', 'param{8}');
%     str = strrep(str, 'param(:,9)', 'param{9}');
%     str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
%     str = strrep(str, 'in5(:,2)', 'ones(ng,1)');
%     if i==7
%         str = '[ng,nc] = size(uh);';
%         fprintf(gid, '%s\n', str);
%         str = ['nch = ' num2str(nch) ';'];
%         fprintf(gid, '%s\n', str);
%         str = ['nd = ' num2str(nd) ';'];
%     end
%     fprintf(gid, '%s\n', str);
%     tline = fgetl(fid);
%     i=i+1;
%     %disp(str)
% end
% str = 'An = reshape(An,ng,nch,nch);';
% fprintf(gid, '%s\n', str);
% str = 'Anm = reshape(Anm,ng,nch,nch,nch);';
% fprintf(gid, '%s\n', str);
% fclose(fid);
% fclose(gid);
% delete('tmp.m');
%
%
% r = uh(1);
% ru = uh(2);
% rv = uh(3);
% rw = uh(4);
% rE = uh(5);
% nx = nl(1);
% ny = nl(2);
% nz = nl(3);
%
% gam1 = gamma-1;
%
% cy = sqrt(nx^2+ny^2);
% sy = nz;
% cz = nx/cy;
% sz = ny/cy;
%
% Tz = [1, 0, 0, 0, 0;...
% 0, cz, sz, 0, 0;...
% 0, -sz, cz, 0, 0;...
% 0, 0, 0, 1, 0;
% 0, 0, 0, 0, 1];
%
% Ty = [1, 0, 0, 0, 0;...
% 0, cy, 0, sy, 0;...
% 0, 0, 1, 0, 0;...
% 0, -sy, 0, cy, 0;
% 0, 0, 0, 0, 1];
%
% T = Ty*Tz;
%
% Tinv = (Tz.')*(Ty.');
%
% U = [r, ru, rv, rw, rE].';
% TU = T*U;
%
% run = TU(2);
% rvn = TU(3);
% rwn = TU(4);
%
% un = run/r;
% vn = rvn/r;
% wn = rwn/r;
%
% P = gam1*(rE - r*(1/2)*(un^2 + vn^2 + wn^2));
% a = sqrt(gamma*P/r);
% H = rE/r + P/r;
%
% A = [ 0 , 1 , 0 , 0, 0 ;...
% gam1*H-un^2-a^2, (3-gamma)*un , -gam1*vn , -gam1*wn, gam1 ;...
% -un*vn , vn , un , 0 , 0 ;...
% -un*wn, wn, 0, un, 0;...
% un*((1/2)*gam1*(un^2+vn^2+wn^2)-H), H-gam1*(un^2), -gam1*un*vn , -gam1*un*wn, gamma*un];
%
% K = [ 1 , 1 , 0 , 0, 1 ;...
% un-a , un , 0 , 0, un+a ;...
% vn , vn , 1 , 0, vn ;...
% wn , wn , 0 , 1, wn ;...
% H - un*a , (1/2)*(un^2+vn^2+wn^2) , vn , wn, H+un*a ];
%
% Kinv = (gam1/(2*a^2))*[ H + (a/gam1)*(un-a) , -(un+a/gam1) , -vn , -wn, 1 ;...
% -2*H + (4/gam1)*a^2 , 2*un , 2*vn , 2*wn, -2 ;...
% -2*(vn*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0, 0 ;...
% -2*wn*a^2/gam1 , 0 , 0 , 2*a^2/gam1, 0
% H - a*(un+a)/gam1 , -un+a/gam1 , -vn , -wn, 1 ];
%
% flux_type = 2;
%
% switch flux_type
% case 2 %%% Lax Friedrich
% An = [ abs(un)+a , 0 , 0 , 0 ;...
% 0 , abs(un)+a , 0 , 0 ;...
% 0 , 0 , abs(un)+a , 0 ;...
% 0 , 0 , 0 , abs(un)+a ];
% case 0 %%% Roe w/o entropy fix
% Lambda = [ un-a , 0 , 0 , 0 , 0 ;...
% 0 , un , 0 , 0 , 0;...
% 0 , 0 , un , 0 , 0;...
% 0 , 0 , 0 , un , 0; ...
% 0 , 0 , 0 , 0 , un+a ];
%
% An = Tinv * (K * Lambda * Kinv) * T;
%
% case 1 %%% Absolute_Roe w/o entropy fix
% Lambda = [ abs(un-a) , 0 , 0 , 0 ;...
% 0 , abs(un) , 0 , 0 ;...
% 0 , 0 , abs(un) , 0 ;...
% 0 , 0 , 0 , abs(un+a) ];
%
% An = K * Lambda * Kinv;
% end
