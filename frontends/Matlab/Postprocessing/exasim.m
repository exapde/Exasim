function [sol,pde,mesh,master,dmd,compilerstr,runstr,res] = exasim(pde,mesh)

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
    % search compilers and set options
    % pde = setcompilers(pde);       
    
    % generate input files and store them in datain folder
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);

    % generate source codes and store them in app folder
    if pde.gencode==1
      %gencode(pde);
      kkgencode(pde);
      compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
    end
        
    runstr = runcode(pde, 1); % run C++ code

    % get solution from output files in dataout folder
    sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');

    % get residual norms from output files in dataout folder
    if pde.saveResNorm
        fileID = fopen('dataout/out_residualnorms0.bin','r'); res = fread(fileID,'double'); fclose(fileID);
        res = reshape(res,4,[])';
    end    
else
    master = cell(nmodels,1);
    dmd = cell(nmodels,1);
    sol = cell(nmodels,1);
    res = cell(nmodels,1);
    
    if pde{1}.gencode==1
      % preprocess and generate code for all PDE models
      for m = 1:nmodels    
          [pde{m},mesh{m},master{m},dmd{m}] = generatecode(pde{m},mesh{m});
      end              
      kkgencodeall(nmodels, pde{m}.exasimpath + "/build/model");
      compilerstr = cmakecompile(pde{1}); % use cmake to compile C++ source codes 
    end
    
    runstr = runcode(pde{1}, nmodels); % run C++ code
    
    % get solution from output files in dataout folder
    for m = 1:nmodels        
        sol{m} = fetchsolution(pde{m},master{m},dmd{m}, pde{m}.buildpath + "/dataout" + num2str(m));        
%         % get residual norms from output files in dataout folder
%         if pde{m}.saveResNorm
%             fileID = fopen(['dataout/out_residualnorms' num2str(m-1) '.bin'],'r'); tm = fread(fileID,'double'); fclose(fileID);
%             res{m} = reshape(tm,4,[])';
%         end            
    end    
end



