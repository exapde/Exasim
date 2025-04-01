syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20
syms x1 x2 x3 av
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10
syms zero one 

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10];
nd=2;
if nd==2    
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12];        
    pg = [x1 x2 x3];       
else
    udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20];        
    pg = [x1 x2 x3];      
end

gam  = param(1);
gam1 = param(1) - 1.0;          
Minf = param(2);
Re   = param(4);
Pr   = param(5);
Tref = param(6);
Twall = param(7);
Re1  = 1/Re;
M2   = Minf^2;
c23  = 2/3;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);

hpar = 1.0;
apar = 10.0;
bpar = 0.5;
cpar = 20.0;

ncu = nd+2;
nc = length(udg);
nch = ncu;

if nd==2                                               
    r    = udg(1);
    ru   = udg(2);
    rv   = udg(3);
    rE   = udg(4);
    rx   = udg(5);
    rux  = udg(6);
    rvx  = udg(7);
    rEx  = udg(8);
    ry   = udg(9);
    ruy  = udg(10);
    rvy  = udg(11);
    rEy  = udg(12);
    av = pg(3);

    r1   = 1./r;
    uv   = ru.*r1;
    vv   = rv.*r1;
    u = uv;
    v = vv;
    E    = rE.*r1;
    q   = 0.5*(uv.*uv+vv.*vv);
    p    = gam1*(rE-r.*q);

    p0 = pinf/2;
    pt = (p-p0).*(atan(3e2*(p-p0))/pi + 1/2) + 1/2 - atan(3e2)/pi + p0;
    dpt = atan(300*p - 300*p0)/pi + (300*(p - p0))./(pi*((300*p - 300*p0).^2 + 1)) + 1/2;
    h    = E+pt.*r1;
    
    fi   = [ru, ru*u+pt, rv*u,   ru*h, ...
            rv, ru*v,   rv*v+pt, rv*h];            
    
  Ts = 110.4;
  T = pt./(gam1*r);  
  Tphy =  (Tref/Tinf) * T;
  muRef = Re1;
  mu = muRef*(Tphy./Tref).^(3/2) .* (Tref + Ts)./(Tphy + Ts);  
  fc = mu*gam/Pr;
  
  ux  = (rux - rx.*uv).*r1;
  vx  = (rvx - rx.*vv).*r1;
  qx  = uv.*ux + vv.*vx;
  px  = gam1*(rEx - rx.*q - r.*qx);
  ptx  = px*dpt;
  Tx  = (1/gam1)*(ptx.*r - pt.*rx).*r1.^2;
     
  uy  = (ruy - ry.*uv).*r1;
  vy  = (rvy - ry.*vv).*r1;
  qy  = uv.*uy + vv.*vy;
  py  = gam1*(rEy - ry.*q - r.*qy);
  pty  = py*dpt;
  Ty  = (1/gam1)*(pty.*r - pt.*ry).*r1.^2;

  txx = c23*mu.*(2*ux - vy);
  txy = mu.*(uy + vx);
  tyy = c23*mu.*(2*vy - ux);

    rHx = rEx+ptx;
    rHy = rEy+pty;           
  
    fv = [0, txx, txy, u*txx + v*txy + fc*Tx,...
          0, txy, tyy, u*txy + v*tyy + fc*Ty];        

    fa = [av.*rx, av.*rux, av.*rvx, av.*(rHx),...
          av.*ry, av.*ruy, av.*rvy, av.*(rHy)];
      
    f1 = fi+fv+fa;  
    
  
%     r1   = 1./r;
%     uv   = ru.*r1;
%     vv   = rv.*r1;
%     E    = rE.*r1;
%     q   = 0.5*(uv.*uv+vv.*vv);
%     p    = gam1*(rE-r.*q);    
% 
%     r1   = 1/udg(1);
%     uv    = udg(2)*r1;
%     vv    = udg(3)*r1;
%     u = uv;
%     v = vv;
%     E    = udg(4)*r1;
%     q    = 0.5*(u*u+v*v);
%     p    = gam1*(udg(4)-udg(1)*q);    
%     
%     pinf = 1/(gam*Minf^2);
%     p0 = pinf/2;
%     pt = (p-p0).*(atan(3e2*(p-p0))/pi + 1/2) + 1/2 - atan(3e2)/pi + p0;
%     dpt = atan(300*p - 300*p0)/pi + (300*(p - p0))./(pi*((300*p - 300*p0).^2 + 1)) + 1/2;
%     h    = E+pt.*r1;
% 
%     fi   = [udg(2), udg(2)*u+pt, udg(3)*u,   udg(2)*h, ...
%             udg(3), udg(2)*v,   udg(3)*v+pt, udg(3)*h];            
% 
%     u_r  = -u.*r1;
%     u_ru =  r1;
%     v_r  = -v.*r1;
%     v_rv =  r1;

%     Ts = 110.4;
%     T = pt./(gam1*udg(1));  
%     Tphy =  (Tref/Tinf) * T;    
%     mu = Re1*(Tphy./Tref).^(3/2) .* (Tref + Ts)./(Tphy + Ts);  
%     fc = mu*gam/Pr;
%     
%     ux  = (udg(6) - udg(5).*uv).*r1;
%     vx  = (udg(7) - udg(5).*vv).*r1;
%     qx  = uv.*ux + vv.*vx;
%     px  = gam1*(udg(8) - udg(5).*q - udg(1).*qx);
%     ptx = dpt*px;
%     Tx  = (1/gam1)*(ptx.*udg(1) - pt.*udg(5)).*r1.^2;
% 
%     uy  = (udg(10) - udg(9).*uv).*r1;
%     vy  = (udg(11) - udg(9).*vv).*r1;
%     qy  = uv.*uy + vv.*vy;
%     py  = gam1*(udg(12) - udg(9).*q - udg(1).*qy);
%     pty = dpt*py;
%     Ty  = (1/gam1)*(pty.*udg(1) - pt.*udg(9)).*r1.^2;
% 
%     txx = c23*mu.*(2*ux - vy);
%     txy = mu.*(uy + vx);
%     tyy = c23*mu.*(2*vy - ux);
%     
%     fv = [0, txx, txy, u*txx + v*txy + fc*Tx,...
%           0, txy, tyy, u*txy + v*tyy + fc*Ty];        
%       
%     rHx = udg(8)+ptx;
%     rHy = udg(12)+pty;           
%     
%     fa = [pg(3)*udg(5), pg(3)*udg(6), pg(3)*udg(7), pg(3)*(rHx),...
%           pg(3)*udg(9), pg(3)*udg(10), pg(3)*udg(11), pg(3)*(rHy)];
        
end


for n=1:1
    if n==1
        f = f1;
        filename1 = ['flux' num2str(nd) 'd2' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['flux' num2str(nd) 'd2' '.m'];
    end
    
    %%% compute Jacobian
    jac_f  = jacobian(f,udg);

    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg);
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg);
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end    
    
    %matlabFunction(f(1),'file','tmp.m','vars',{pg,udg,param,time,[zero one]},'outputs', {'f'});
    
    % generate a temporary matlab file
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{pg,udg,param,time,[zero one]},'outputs', {'f','f_udg'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename1,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
        str = strrep(str, 'in1', 'pg');
        str = strrep(str, 'IN1', 'PG');        
        str = strrep(str, 'in2', 'udg');
        str = strrep(str, 'IN2', 'UDG');        
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
        str = strrep(str, 'param(:,10)', 'param{10}');     
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
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

    str = 'f = reshape(f,ng,nch,nd);';
    fprintf(gid, '%s\n', str);                  
    str = 'f_udg = reshape(f_udg,ng,nch,nd,nc);';
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    %delete('tmp.m');
end

