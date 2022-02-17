

filein = string(filepath) + '/' + string(filename);                 
run(filein);
if iscell(sol)
    for j = 1:length(sol)
        filesol = string(filepath) + '/' + 'sol' + num2str(j) + '.bin'; 
        fileID = fopen(filesol,'w');
        fwrite(fileID,sol{j}(:),'double','native');
        fclose(fileID);      
        rmdir(string(filepath) + '/' + 'dataout' + num2str(j), 's');         
        rmdir(string(filepath) + '/' + 'datain'  + num2str(j), 's');         
    end            
    rmdir(string(filepath) + '/' + 'app', 's');                            
else
    filesol = string(filepath) + '/' + 'sol.bin'; 
    fileID = fopen(filesol,'w');
    fwrite(fileID,sol(:),'double','native');
    fclose(fileID);        
    rmdir(string(filepath) + '/' + 'dataout', 's');         
    rmdir(string(filepath) + '/' + 'datain', 's');         
    rmdir(string(filepath) + '/' + 'app', 's');                            
end

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
cmdWindow = desktop.getClient('Command Window');
cmdWindowScrollPane = cmdWindow.getComponent(0);
cmdWindowScrollPaneViewport = cmdWindowScrollPane.getComponent(0);
cmdTextUI = cmdWindowScrollPaneViewport.getComponent(0);
cmdText = cmdTextUI.getText;

fileID = fopen(string(filepath) + '/' + 'log.txt','w');
fprintf(fileID, '%c', char(cmdText));        
fclose(fileID);                
clc;        


