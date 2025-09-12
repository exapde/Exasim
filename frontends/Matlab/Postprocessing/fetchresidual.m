function res = fetchresidual(pde)

if iscell(pde)
    nmodels = length(pde);
    if nmodels==1
        error("Number of PDE models must be greater than 1");
    end
else
    nmodels = 1;
end

res = [];
if nmodels==1    
    % get residual norms from output files in dataout folder
    if pde.saveResNorm
        fileID = fopen('dataout/out_residualnorms0.bin','r'); res = fread(fileID,'double'); fclose(fileID);
        res = reshape(res,4,[])';
    end    
else
    res = cell(nmodels,1);    
    % get solution from output files in dataout folder
    for m = 1:nmodels                
        % get residual norms from output files in dataout folder
        if pde{m}.saveResNorm
            fileID = fopen(['dataout/out_residualnorms' num2str(m-1) '.bin'],'r'); tm = fread(fileID,'double'); fclose(fileID);
            res{m} = reshape(tm,4,[])';
        end            
    end    
end



