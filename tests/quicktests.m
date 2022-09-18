cdir = pwd;
ii = strfind(cdir, "Exasim"); 
testdir = cdir(1:(ii+5)) + "/tests";
cd(testdir);
mfiles = dir('**/*.m');
n = length(mfiles);
success = ones(n,1);

clc;
for i = 1:n
    if string(mfiles(i).name) == "pdeapp.m"    
        s = string(mfiles(i).folder);
        if contains(s,"Cone")==1 || contains(s,"hypersoniccylinder")==1 || contains(s,"TaylorGreenVortex")==1 || contains(s,"naca0012")==1 || contains(s,"Scattering")==1 || contains(s,"Ring2D")==1
            disp(s);
        else
        filein = string(mfiles(i).folder) + '/' + string(mfiles(i).name);                 
        run(filein);                
        if iscell(sol)
            for j = 1:length(sol)
                filename = string(mfiles(i).folder) + '/' + 'sol' + num2str(j) + '.bin'; 
                fileID = fopen(filename,'r');                
                tmp = fread(fileID,'double');  
                fclose(fileID);      
                rmdir(string(mfiles(i).folder) + '/' + 'dataout' + num2str(j), 's');         
                rmdir(string(mfiles(i).folder) + '/' + 'datain'  + num2str(j), 's');         
                err = max(abs(tmp(:)-sol{j}(:)));
                if err<=1e-10
                    disp("solutions match");
                else
                    warning("solutions do not match");
                end
            end            
            rmdir(string(mfiles(i).folder) + '/' + 'app', 's');                            
        else
            filename = string(mfiles(i).folder) + '/' + 'sol.bin'; 
            fileID = fopen(filename,'r');
            tmp = fread(fileID,'double');                        
            fclose(fileID);        
            rmdir(string(mfiles(i).folder) + '/' + 'dataout', 's');         
            rmdir(string(mfiles(i).folder) + '/' + 'datain', 's');         
            rmdir(string(mfiles(i).folder) + '/' + 'app', 's');       
            if contains(string(mfiles(i).folder),"NavierStokes")
                if size(mesh.p,1)==2
                    sol = sol(:,1:4,:,end);
                elseif size(mesh.p,1)==3
                    sol = sol(:,1:5,:,end);
                end                
            end
            if contains(string(mfiles(i).folder),"WaveEquation")                
                sol = sol(:,:,:,end);                
            end
            err = max(abs(tmp(:)-sol(:)));
            if err<=1e-10
                disp("solutions match");
            else
                disp(filein);
                warning("solutions do not match. Press any key to continue");                
                success(i) = 0;
                pause;                
            end
        end
        
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        cmdWindow = desktop.getClient('Command Window');
        cmdWindowScrollPane = cmdWindow.getComponent(0);
        cmdWindowScrollPaneViewport = cmdWindowScrollPane.getComponent(0);
        cmdTextUI = cmdWindowScrollPaneViewport.getComponent(0);
        cmdText = cmdTextUI.getText;        
        q = char(cmdText);
        ii = strfind(q,'PTC Iteration:'); ii = ii(1);
        q = q(ii:end);
        
        fileID = fopen(string(mfiles(i).folder) + '/' + 'log.txt','r');
        s = fscanf(fileID,'%c');        
        fclose(fileID);                        
        ii = strfind(s,'PTC Iteration:'); ii = ii(1);
        s = s(ii:end);
        
        if string(q(1:length(s))) == string(s)
            disp("log files match");
        else
            if success(i) == 0
                success(i) = -2;
            else
                success(i) = -1;
            end
            disp(filein);
            warning("log files do not match. Press any key to continue");            
            pause;
        end        
        clc;        
        end
    end    
end

for i = 1:n      
    if string(mfiles(i).name) == "pdeapp.m"        
        if success(i) == 0
            disp(filein);
            warning("solutions do not match");                            
        elseif success(i) == -1
            disp(filein);
            warning("log files do not match");                        
        elseif success(i) == -2
            disp(filein);
            warning("solutions and log files do not match");                                    
        end
    end
end


