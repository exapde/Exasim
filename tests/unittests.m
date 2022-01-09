cdir = pwd;
ii = strfind(cdir, "Exasim"); 
testdir = cdir(1:(ii+5)) + "/tests";
cd(testdir);
mfiles = dir('**/*.m');
n = length(mfiles);

clc;
for i = 1:n
    if string(mfiles(i).name) == "pdeapp.m"        
        filein = string(mfiles(i).folder) + '/' + string(mfiles(i).name);                 
        run(filein);                
        if iscell(sol)
            for j = 1:length(sol)
                filename = string(mfiles(i).folder) + '/' + 'sol' + num2str(j) + '.bin'; 
                fileID = fopen(filename,'w');
                fwrite(fileID,sol{j}(:),'double','native');
                fclose(fileID);      
                rmdir(string(mfiles(i).folder) + '/' + 'dataout' + num2str(j), 's');         
                rmdir(string(mfiles(i).folder) + '/' + 'datain'  + num2str(j), 's');         
            end            
            rmdir(string(mfiles(i).folder) + '/' + 'app', 's');                            
        else
            filename = string(mfiles(i).folder) + '/' + 'sol.bin'; 
            fileID = fopen(filename,'w');
            fwrite(fileID,sol(:),'double','native');
            fclose(fileID);        
            rmdir(string(mfiles(i).folder) + '/' + 'dataout', 's');         
            rmdir(string(mfiles(i).folder) + '/' + 'datain', 's');         
            rmdir(string(mfiles(i).folder) + '/' + 'app', 's');                            
        end
        
        desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
        cmdWindow = desktop.getClient('Command Window');
        cmdWindowScrollPane = cmdWindow.getComponent(0);
        cmdWindowScrollPaneViewport = cmdWindowScrollPane.getComponent(0);
        cmdTextUI = cmdWindowScrollPaneViewport.getComponent(0);
        cmdText = cmdTextUI.getText;
        
        fileID = fopen(string(mfiles(i).folder) + '/' + 'log.txt','w');
        fprintf(fileID, '%c', char(cmdText));        
        fclose(fileID);                
        clc;        
    end    
end


