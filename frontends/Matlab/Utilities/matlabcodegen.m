
syms zero one 

% generate matlab code for the flux function
if exist('f','var') && exist('udg','var')
    
    %%% compute Jacobian
    jac_f  = jacobian(f,udg);

    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg,'legacy');
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg,'legacy');
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

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
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', join(str));                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', join(str));                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', join(str));                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    str = 'f = reshape(f,ng,nch,nd);';
    fprintf(gid, '%s\n', join(str));                  
    str = 'f_udg = reshape(f_udg,ng,nch,nd,nc);';
    fprintf(gid, '%s\n', join(str));                  
    
    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end


% generate matlab code for the source function
if exist('s','var') && exist('udg','var')
    jac_s  = jacobian(s,udg);

    for ii = 1:size(s,1)
        for jj = 1:size(s,2)
            temp = ismember(symvar(s(ii,jj)),udg,'legacy');
            if s(ii,jj)==0, s(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, s(ii,jj) = (s(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_s,1)
        for jj = 1:size(jac_s,2)
            temp = ismember(symvar(jac_s(ii,jj)),udg,'legacy');
            if jac_s(ii,jj)==0, jac_s(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_s(ii,jj) = (jac_s(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    % generate a temporary matlab file
    matlabFunction(s(:),jac_s(:),'file','tmp.m',...
        'vars',{pg,udg,param,time,[zero one]},'outputs', {'s','s_udg'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename2,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename2,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename2,'.m','')));
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
        str = strrep(str, 'in5(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in5(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', join(str));                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', join(str));                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', join(str));                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    str = 's = reshape(s,ng,nch);';
    fprintf(gid, '%s\n', join(str));                  
    str = 's_udg = reshape(s_udg,ng,nch,nc);';
    fprintf(gid, '%s\n', join(str));                  
    
    fclose(fid);
    fclose(gid);
    %delete('tmp.m');
elseif exist('filename2','var')
    gid = fopen(filename2,'wt');
    
    str = ['function [s,s_udg] = ' strrep(filename2,'.m','') '(pg,udg,param,time)'];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%' upper(filename2)];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%    [s,s_udg] = ' strrep(filename2,'.m','') '(pg,udg,param,time)'];
    fprintf(gid, '%s\n', upper(join(str)));                  
    str = '[ng,nc] = size(udg);';
    fprintf(gid, '%s\n', join(str));                  
    str = ['nch = ' num2str(nch) ';'];
    fprintf(gid, '%s\n', join(str));                          
    str = 's = zeros(ng,nch);';
    fprintf(gid, '%s\n', join(str));                  
    str = 's_udg = zeros(ng,nch,nc);';
    fprintf(gid, '%s\n', join(str));                  
    
    fclose(gid);    
end

% generate matlab code for the fhat function
if exist('fh','var') && exist('udg','var') && exist('uh','var')
    
    jac_fh   = jacobian(fh,udg);
    jac_fhuh = jacobian(fh,uh);

    for ii = 1:size(fh,1)
        for jj = 1:size(fh,2)
            temp = ismember(symvar(fh(ii,jj)),udg,'legacy');
            if fh(ii,jj)==0, fh(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, fh(ii,jj) = (fh(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_fh,1)
        for jj = 1:size(jac_fh,2)
            temp = ismember(symvar(jac_fh(ii,jj)),udg,'legacy');
            if jac_fh(ii,jj)==0, jac_fh(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_fh(ii,jj) = (jac_fh(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_fhuh,1)
        for jj = 1:size(jac_fhuh,2)
            temp = ismember(symvar(jac_fhuh(ii,jj)),udg,'legacy');
            if jac_fhuh(ii,jj)==0, jac_fhuh(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_fhuh(ii,jj) = (jac_fhuh(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end    
    
    % generate a temporary matlab file    
    matlabFunction(fh(:),jac_fh(:),jac_fhuh(:),'file','tmp.m',...
        'vars',{nl,pg,udg,uh,param,time,[zero one]},'outputs', {'fh','fh_udg','fh_uh'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename3,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename3,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename3,'.m','')));
        str = strrep(str, 'in1', 'nl');
        str = strrep(str, 'IN1', 'NL');        
        str = strrep(str, 'in2', 'pg');
        str = strrep(str, 'IN2', 'PG');        
        str = strrep(str, 'in3', 'udg');
        str = strrep(str, 'IN3', 'UDG');        
        str = strrep(str, 'in4', 'uh');
        str = strrep(str, 'IN4', 'UH');                
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
        str = strrep(str, 'in7(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in7(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', join(str));                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', join(str));                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', join(str));                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    str = 'fh = reshape(fh,ng,nch);';
    fprintf(gid, '%s\n', join(str));                  
    str = 'fh_udg = reshape(fh_udg,ng,nch,nc);';
    fprintf(gid, '%s\n', join(str));                  
    str = 'fh_uh = reshape(fh_uh,ng,nch,nch);';
    fprintf(gid, '%s\n', join(str));                  
    
    fclose(fid);
    fclose(gid);
    delete('tmp.m');
elseif exist('filename3','var')
    gid = fopen(filename3,'wt');
    
    str = ['function [fh,fh_udg,fh_uh] = ' strrep(filename3,'.m','') '(nl,pg,udg,uh,param,time)'];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%' upper(filename3)];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%    [fh,fh_udg,fh_uh] = ' strrep(filename3,'.m','') '(nl,pg,udg,uh,param,time)'];
    fprintf(gid, '%s\n', upper(join(str)));                  
    str = '[ng,nc] = size(udg);';
    fprintf(gid, '%s\n', join(str));                  
    str = 'nch = size(uh,2);';
    fprintf(gid, '%s\n', join(str));                          
    str = 'tau = param{end};';
    fprintf(gid, '%s\n', join(str));                          
    
    str = 'u = udg(:,1:nch);';
    fprintf(gid, '%s\n', join(str));                             
    str = '[f,f_uhq] = flux(pg,[uh,udg(:,nch+1:end)],param,time);';
    fprintf(gid, '%s\n', join(str));                              
    
    str = 'fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + tau*(u-uh),[2 1]);';
    fprintf(gid, '%s\n', join(str));                          
    str = 'fh_u = tau*ones(ng,nch);';
    fprintf(gid, '%s\n', join(str));                          
    
    str = 'fh_uhq = mapContractK(f_uhq,nl,[2 4],3,1,2,[],1);';
    fprintf(gid, '%s\n', join(str));                          
    str = 'fh_q = permute(fh_uhq(:,nch+1:nc,:),[3 1 2]);';
    fprintf(gid, '%s\n', join(str));                          
    str = 'fh_udg = cat(3,fh_u,fh_q);';
    fprintf(gid, '%s\n', join(str));         
    
    str = 'fh_uh = permute(fh_uhq(:,1:nch,:),[3 1 2]) - tau;';    
    fprintf(gid, '%s\n', join(str));         
        
    fclose(gid);        
end

% generate matlab code for the fbou function
if exist('fb','var') && exist('udg','var') && exist('uh','var')
    
    gid = fopen(filename4,'wt');
    
    str = ['function [fb,fb_udg,fb_uh] = ' strrep(filename4,'.m','') '(ib,uinf,nl,pg,udg,uh,param,time)'];        
    fprintf(gid, '%s\n', join(str));                  
    str = ['%' upper(filename4)];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%    [fb,fb_udg,fb_uh] = ' strrep(filename4,'.m','') '(ib,uinf,nl,pg,udg,uh,param,time)'];
    fprintf(gid, '%s\n', upper(join(str)));                  
    str = '[ng,nc] = size(udg);';
    fprintf(gid, '%s\n', join(str));                  
    str = ['nch = ' num2str(nch) ';'];
    fprintf(gid, '%s\n', join(str));                  
    str = ['nd = ' num2str(nd) ';'];
    fprintf(gid, '%s\n', join(str));      
    str = 'switch (ib)';
    fprintf(gid, '%s\n', join(str));  
    
    for k=1:length(fb)
        str = ['    case ' num2str(k)];
        fprintf(gid, '%s\n', join(str));  
        if strcmp(fb{k},'dirichlet')==1
            str = '        fb = uinf-uh;';
            fprintf(gid, '%s\n', join(str));  
            str = '        fb_udg = zeros(ng,nch,nc);';
            fprintf(gid, '%s\n', join(str));  
            str = '        fb_uh = -ones(ng,nch,nch);';
            fprintf(gid, '%s\n', join(str));  
        elseif strcmp(fb{k},'neumann')==1                           
            str = '        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);';
            fprintf(gid, '%s\n', join(str));  
            str = '        fb = fb - uinf;';
            fprintf(gid, '%s\n', join(str));  
        else
            fbk = fb{k};
            jac_fh   = jacobian(fbk,udg);
            jac_fhuh = jacobian(fbk,uh);

            for ii = 1:size(fb,1)
                for jj = 1:size(fbk,2)
                    temp = ismember(symvar(fbk(ii,jj)),udg,'legacy');
                    if fbk(ii,jj)==0, fbk(ii,jj) = zero; %%% No dependency on anything
                    elseif isempty(temp) || sum(temp)==0, fbk(ii,jj) = (fbk(ii,jj))*one; %%% No dependency on state vars
                    end
                end
            end

            for ii = 1:size(jac_fh,1)
                for jj = 1:size(jac_fh,2)
                    temp = ismember(symvar(jac_fh(ii,jj)),udg,'legacy');
                    if jac_fh(ii,jj)==0, jac_fh(ii,jj) = zero; %%% No dependency on anything
                    elseif isempty(temp) || sum(temp)==0, jac_fh(ii,jj) = (jac_fh(ii,jj))*one; %%% No dependency on state vars
                    end
                end
            end

            for ii = 1:size(jac_fhuh,1)
                for jj = 1:size(jac_fhuh,2)
                    temp = ismember(symvar(jac_fhuh(ii,jj)),udg,'legacy');
                    if jac_fhuh(ii,jj)==0, jac_fhuh(ii,jj) = zero; %%% No dependency on anything
                    elseif isempty(temp) || sum(temp)==0, jac_fhuh(ii,jj) = (jac_fhuh(ii,jj))*one; %%% No dependency on state vars
                    end
                end
            end    

            % generate a temporary matlab file    
            matlabFunction(fbk(:),jac_fh(:),jac_fhuh(:),'file','tmp.m','vars',{uinf,nl,pg,udg,uh,param,time,[zero one]},'outputs', {'fb','fb_udg','fb_uh'});
                        
            % open the file and modify it
            fid = fopen('tmp.m','r');    

            tline = fgetl(fid); 
            i=1;       
            while ischar(tline)        
                if i>7
                str = strrep(tline, 'in1', 'uinf');
                str = strrep(str, 'IN1', 'uinf');            
                str = strrep(str, 'in2', 'nl');
                str = strrep(str, 'IN2', 'NL');        
                str = strrep(str, 'in3', 'pg');
                str = strrep(str, 'IN3', 'PG');        
                str = strrep(str, 'in4', 'udg');
                str = strrep(str, 'IN4', 'UDG');        
                str = strrep(str, 'in5', 'uh');
                str = strrep(str, 'IN5', 'UH');                
                str = strrep(str, 'in6', 'param');
                str = strrep(str, 'IN6', 'PARAM');                                
                str = strrep(str, 'param(:,1)', 'param{1}'); 
                str = strrep(str, 'param(:,2)', 'param{2}'); 
                str = strrep(str, 'param(:,3)', 'param{3}'); 
                str = strrep(str, 'param(:,4)', 'param{4}'); 
                str = strrep(str, 'param(:,5)', 'param{5}'); 
                str = strrep(str, 'param(:,6)', 'param{6}'); 
                str = strrep(str, 'param(:,7)', 'param{7}'); 
                str = strrep(str, 'param(:,8)', 'param{8}'); 
                str = strrep(str, 'param(:,9)', 'param{9}');         
                str = strrep(str, 'in8(:,1)', 'zeros(ng,1)');
                str = strrep(str, 'in8(:,2)', 'ones(ng,1)');
                str = ['        ' str];
                fprintf(gid, '%s\n', join(str));                                  
                end
                tline = fgetl(fid);
                i=i+1;                                
            end

            str = '        fb = reshape(fb,ng,nch);';
            fprintf(gid, '%s\n', join(str));                  
            str = '        fb_udg = reshape(fb_udg,ng,nch,nc);';
            fprintf(gid, '%s\n', join(str));                  
            str = '        fb_uh = reshape(fb_uh,ng,nch,nch);';
            fprintf(gid, '%s\n', join(str));                  

            fclose(fid);
            %delete('tmp.m');            
        end    
    end
    
    str = '    otherwise';
    fprintf(gid, '%s\n', join(str));    
    str = '         error(''unknown boundary type'');';
    fprintf(gid, '%s\n', join(str));    
    str = 'end';
    fprintf(gid, '%s\n', join(str));  
    fclose(gid);  
else
    gid = fopen(filename4,'wt');
    
    str = ['function [fb,fb_udg,fb_uh] = ' strrep(filename4,'.m','') '(ib,uinf,nl,pg,udg,uh,param,time)'];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%' upper(filename4)];
    fprintf(gid, '%s\n', join(str));                  
    str = ['%    [fb,fb_udg,fb_uh] = ' strrep(filename4,'.m','') '(ib,uinf,nl,pg,udg,uh,param,time)'];
    fprintf(gid, '%s\n', upper(join(str)));                  
    str = '[ng,nc] = size(udg);';
    fprintf(gid, '%s\n', join(str));                  
    str = ['nch = ' num2str(nch) ';'];
    fprintf(gid, '%s\n', join(str));                      
    str = 'switch (ib)';
    fprintf(gid, '%s\n', join(str));  
        
    str = ['    case ' num2str(1)];
    fprintf(gid, '%s\n', join(str));          
    str = '        fb = uinf-uh;';
    fprintf(gid, '%s\n', join(str));  
    str = '        fb_udg = zeros(ng,nch,nc);';
    fprintf(gid, '%s\n', join(str));  
    str = '        fb_uh = -ones(ng,nch,nch);';
    fprintf(gid, '%s\n', join(str));  
        
    str = ['    case ' num2str(2)];    
    fprintf(gid, '%s\n', join(str));          
    str = '        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);';
    fprintf(gid, '%s\n', join(str));  
    str = '        fb = fb - uinf;';
    fprintf(gid, '%s\n', join(str));  
    
    str = '    otherwise';
    fprintf(gid, '%s\n', join(str));    
    str = '         error(''unknown boundary type'');';
    fprintf(gid, '%s\n', join(str));    
    str = 'end';
    fprintf(gid, '%s\n', join(str));  
    fclose(gid);      
end


% % generate matlab code for the stabilization function
% if exist('st','var') && exist('uh','var')
%         
%     jac_stuh = jacobian(st,uh);
% 
%     for ii = 1:size(st,1)
%         for jj = 1:size(st,2)
%             temp = ismember(symvar(st(ii,jj)),udg);
%             if st(ii,jj)==0, st(ii,jj) = zero; %%% No dependency on anything
%             elseif isempty(temp) || sum(temp)==0, st(ii,jj) = (st(ii,jj))*one; %%% No dependency on state vars
%             end
%         end
%     end
% 
%     for ii = 1:size(jac_stuh,1)
%         for jj = 1:size(jac_stuh,2)
%             temp = ismember(symvar(jac_stuh(ii,jj)),udg);
%             if jac_stuh(ii,jj)==0, jac_stuh(ii,jj) = zero; %%% No dependency on anything
%             elseif isempty(temp) || sum(temp)==0, jac_stuh(ii,jj) = (jac_stuh(ii,jj))*one; %%% No dependency on state vars
%             end
%         end
%     end    
%     
%     % generate a temporary matlab file    
%     matlabFunction(st(:),jac_stuh(:),'file','tmp.m',...
%         'vars',{nl,pg,uh,param,time,[zero one]},'outputs', {'st','st_uh'});
% 
%     % open the file and modify it
%     fid = fopen('tmp.m','r');
%     gid = fopen(filename5,'wt');
% 
%     tline = fgetl(fid); 
%     i=1;       
%     while ischar(tline)        
%         str = strrep(tline, 'tmp', strrep(filename5,'.m',''));
%         str = strrep(str, 'TMP', upper(strrep(filename5,'.m','')));
%         str = strrep(str, 'in1', 'nl');
%         str = strrep(str, 'IN1', 'NL');                        
%         str = strrep(str, 'in2', 'pg');
%         str = strrep(str, 'IN2', 'PG');                
%         str = strrep(str, 'in3', 'uh');
%         str = strrep(str, 'IN3', 'UH');                
%         str = strrep(str, 'in4', 'param');
%         str = strrep(str, 'IN4', 'PARAM');                
%         str = strrep(str, ',in6)', ')');                
%         str = strrep(str, ',IN6)', ')');    
%         str = strrep(str, 'param(:,1)', 'param{1}'); 
%         str = strrep(str, 'param(:,2)', 'param{2}'); 
%         str = strrep(str, 'param(:,3)', 'param{3}'); 
%         str = strrep(str, 'param(:,4)', 'param{4}'); 
%         str = strrep(str, 'param(:,5)', 'param{5}'); 
%         str = strrep(str, 'param(:,6)', 'param{6}'); 
%         str = strrep(str, 'param(:,7)', 'param{7}'); 
%         str = strrep(str, 'param(:,8)', 'param{8}'); 
%         str = strrep(str, 'param(:,9)', 'param{9}');         
%         str = strrep(str, 'in6(:,1)', 'zeros(ng,1)');
%         str = strrep(str, 'in6(:,2)', 'ones(ng,1)');        
%         if i==7
%             str = '[ng,nc] = size(udg);';
%             fprintf(gid, '%s\n', join(str));                  
%             str = ['nch = ' num2str(nch) ';'];
%             fprintf(gid, '%s\n', join(str));                  
%             str = ['nd = ' num2str(nd) ';'];
%         end
%         fprintf(gid, '%s\n', join(str));                  
%         tline = fgetl(fid);        
%         i=i+1;
%         %disp(str)
%     end
% 
%     str = 'st = reshape(st,ng,nch,nch);';
%     fprintf(gid, '%s\n', join(str));                      
%     str = 'st_uh = reshape(st_uh,ng,nch,nch,nch);';
%     fprintf(gid, '%s\n', join(str));                  
%     
%     fclose(fid);
%     fclose(gid);
%     delete('tmp.m');
% end
% 

