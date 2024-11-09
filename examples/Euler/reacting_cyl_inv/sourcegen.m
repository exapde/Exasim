% syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24
% syms w1 w2
% syms dw1_du1 dw1_du2 dw1_du3 dw1_du4 dw1_du5 dw1_du6 dw1_du7 dw1_du8
syms x1 x2 x3 av
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13
syms zero one 

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13];
nd=2;

% udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24]; 
udg = sym("u",[1 24]); % (ns+nd+1)*(nd+1)
wdg = sym("w",[1 2]); % needs to be at least 2
% dw_du = sym("dw_du",[2 8]);  
% wdg = [w1 w2];
% dw_du = [dw1_du1 dw1_du2 dw1_du3 dw1_du4 dw1_du5 dw1_du6 dw1_du7 dw1_du8];
pg = [x1 x2 x3];       

ns = 5;
ncu = ns + nd + 1;
nc = length(udg);
nch = ncu;

% Nondimensional params
rho_scale   = param(1);
u_scale     = param(2);
rhoe_scale  = param(3);
T_scale     = param(4);
mu_scale    = param(5);
kappa_scale = param(6);
cp_scale    = param(7);
L_scale     = param(8);
% Ec          = param(9);
Ec = 1.0;
Pr          = param(10);
Re          = param(11);
T_wall = param(12);
omega_scale = rho_scale * u_scale / L_scale;


i_T = 1;
i_P = 2;

kinetics_params = kinetics();
[species_thermo_structs, Mw, RU] = thermodynamicsModels();
Mw = Mw';

% Mutation outputs
rho_i_dim = udg(1:ns) * rho_scale;
T = wdg(i_T) * T_scale;
omega_i = netProductionRatesTotal(rho_i_dim, T, Mw, kinetics_params, species_thermo_structs);

f = sym(zeros(8,1));

f(1:ns) = omega_i / omega_scale;

for n=1
    % filename1 = ['densityWallTemp2Energy' num2str(nd) 'd' '.m'];
    filename1 = 'source2d_gen.m';
    %%% compute Jacobian
    jac_f = jacobian(f,udg);
    jac_w = jacobian(f,wdg);
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
    jac_f(:,1:ncu) = jac_f(:,1:ncu) + jac_w * dw_du;
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{pg,udg,wdg,dw_du,param,time,[zero one]},'outputs', {'f','f_udg'});

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
        str = strrep(str, 'in3', 'wdg');
        str = strrep(str, 'IN3', 'WDG');     
        str = strrep(str, 'in4', 'dw_du');
        str = strrep(str, 'IN4', 'DW_DU');
        str = strrep(str, 'in5', 'param');
        str = strrep(str, 'IN5', 'PARAM');     
        str = strrep(str, ',in7)', ')');                
        str = strrep(str, ',IN7)', ')');    
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
        str = strrep(str, 'param(:,11)', 'param{11}');     
        str = strrep(str, 'in7(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in7(:,2)', 'ones(ng,1)');        
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

    str = 'f = reshape(f,ng,nch);';
    fprintf(gid, '%s\n', str);                  
    str = 'f_udg = reshape(f_udg,ng,nch,nc);';
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end


