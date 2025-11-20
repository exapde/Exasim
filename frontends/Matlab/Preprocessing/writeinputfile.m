function writeinputfile(filename, pde, mesh)

fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file %s for writing.', filename);
end

pde.datapath = pwd();
if pde.hybrid==0
  pde.discretization = "ldg";
elseif pde.hybrid==1
  pde.discretization = "hdg";  
end
pde.NewtonIter =  pde.NLiter;
pde.NewtonTol =  pde.NLtol;
pde.GMRESiter =  pde.linearsolveriter;
pde.GMREStol =  pde.linearsolvertol;
pde.dt = pde.dt(:)';
pde.tau = pde.tau(:)';
pde.physicsparam = pde.physicsparam(:)';
pde.externalparam = pde.externalparam(:)';

pde.boundaryconditions = mesh.boundarycondition(:)';
pde.boundaryexpressions = convertHandlesToStrings(mesh.boundaryexpr);
pde.curvedboundaries = mesh.curvedboundary;
pde.curvedboundaryexprs = mesh.curvedboundaryexpr;

if isfield(pde, 'curvedboundaries') == 0 || isempty(pde.curvedboundaries)
  pde.curvedboundaries = 0*pde.boundaryconditions;
  pde.curvedboundaryexprs = repmat("", [1 length(pde.boundaryconditions)]);
else
  pde.curvedboundaryexprs = convertHandlesToStrings(pde.curvedboundaryexprs);
end

% mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:)};
% mesh.periodicexpr = {2, @(p) p([2 3],:), 4, @(p) p([2 3],:); 1, @(p) p([1 3],:), 3, @(p) p([1 3],:); 5, @(p) p([1 2],:), 6, @(p) p([1 2],:)};

if isfield(mesh, 'periodicboundary') == 0 || isempty(mesh.periodicboundary)
  pde.periodicboundaries1 = [];
  pde.periodicboundaries2 = [];
  pde.periodicexprs1 = [];
  pde.periodicexprs2 = [];
else  
  periodicexpr = mesh.periodicboundary;
  m = size(periodicexpr,1);
  pde.periodicboundaries1 = zeros(1, m);
  pde.periodicboundaries2 = zeros(1, m);
  pde.periodicexprs1 = strings(1,m);
  pde.periodicexprs2 = strings(1,m);
  for j = 1:m
    pde.periodicboundaries1(j) = periodicexpr{j,1};
    pde.periodicboundaries1(j) = periodicexpr{j,3};    
    tm = convertHandlesToStrings({periodicexpr{j,2}});    
    if tm == "xy"
      pde.periodicexprs1(2*j-1:2*j) = ["x", "y"];
    elseif tm == "xz"
      pde.periodicexprs1(2*j-1:2*j) = ["x", "z"];  
    elseif tm == "yz"
      pde.periodicexprs1(2*j-1:2*j) = ["y", "z"];     
    else
      pde.periodicexprs1(j) = tm;
    end    
    tm = convertHandlesToStrings({periodicexpr{j,4}});    
    if tm == "xy"
      pde.periodicexprs2(2*j-1:2*j) = ["x", "y"];
    elseif tm == "xz"
      pde.periodicexprs2(2*j-1:2*j) = ["x", "z"];  
    elseif tm == "yz"
      pde.periodicexprs2(2*j-1:2*j) = ["y", "z"];     
    else
      pde.periodicexprs2(j) = tm;
    end    
  end
end
  
if isfield(mesh, 'interfacecondition')
  pde.interfaceconditions = mesh.interfacecondition;
else
  pde.interfaceconditions = [];
end

if isfield(pde, 'meshfile') == 0
  pde.meshfile = "grid.bin";  
  writebin(pde.meshfile, [size(mesh.p) size(mesh.t) mesh.p(:)' mesh.t(:)']);
end

if isfield(mesh, 'dgnodes') == 1
  pde.xdgfile = "xdg.bin";  
  writebin(pde.xdgfile, [size(mesh.dgnodes) mesh.dgnodes(:)']);
end

if isfield(mesh, 'udg') == 1
  pde.udgfile = "udg.bin";  
  writebin(pde.udgfile, [size(mesh.udg) mesh.udg(:)']);
end

if isfield(mesh, 'vdg') == 1
  pde.vdgfile = "vdg.bin";  
  writebin(pde.vdgfile, [size(mesh.vdg) mesh.vdg(:)']);
end

if isfield(mesh, 'wdg') == 1
  pde.wdgfile = "wdg.bin";  
  writebin(pde.wdgfile, [size(mesh.wdg) mesh.wdg(:)']);
end

%fprintf(fid, 'Input file saved on: %s\n\n', datestr(now));
fields = fieldnames(pde);
requiredKeys = ["exasimpath", "datapath", "model", "modelfile", "meshfile", "xdgfile",...
  "udgfile", "vdgfile", "wdgfile", "discretization",...
  "platform", "mpiprocs", "debugmode", "runmode", "modelnumber", "porder", "pgauss","torder",...
  "nstage","ncu", "ncw", "neb", "nfb", "NewtonIter", "NewtonTol", "GMRESiter", "GMRESrestart",...
  "GMREStol", "GMRESortho","ppdegree","RBdim", "matvecorder", "matvectol", "precMatrixType",...
  "preconditioner", "time", "tau", "dt", "physicsparam", "externalparam", "boundaryconditions",...
  "boundaryexpressions", "curvedboundaries", "curvedboundaryexprs", "periodicboundaries1",...
  "periodicexprs1", "periodicboundaries2","periodicexprs2", "interfaceconditions"];
  
for i = 1:length(requiredKeys)
    key = requiredKeys(i);    
    
    match = 0;
    for j = 1:length(fields)
      if (key == fields(j)) 
        match = 1;
      end
    end
      
    if match == 1
       value = pde.(key);
      if isnumeric(value) || islogical(value)        
        if isempty(value)
          fprintf(fid, '%s = [];\n', key);       
        elseif length(value)==1
          if key == "dt" || key == "tau" || key == "physicsparam" || key == "externalparam" || key == "boundaryconditions" || key == "curvedboundaries" || key == "periodicboundaries1" || key == "periodicboundaries2"
            fprintf(fid, '%s = [%s];\n', key, mat2str(value));       
          else
            fprintf(fid, '%s = %s;\n', key, mat2str(value));       
          end
        else                       
          tm = num2str(value(1));
          for k = 2:length(value)
            tm = tm +  ", " + num2str(value(k));            
          end          
          tm = "[" + tm + "]";          
          fprintf(fid, '%s = %s;\n', key, char(tm));         
        end
      elseif ischar(value)
          if ~((key == "exasimpath") || (key == "datapath"))           
            fprintf(fid, '%s = \"%s\";\n', key, value);  
          end
      elseif isstring(value)          
        if length(value)==1
          if key == "boundaryexpressions" || key == "curvedboundaryexprs" || key == "periodicexprs1" || key == "periodicexprs2"
            fprintf(fid, '%s = [\"%s\"];\n', key, char(value));     
          elseif (key == "modelfile")  
            fprintf(fid, '%s = \"%s\";\n', key, char(value + ".txt"));     
          elseif (key == "exasimpath") || (key == "datapath")          
          else
            fprintf(fid, '%s = \"%s\";\n', key, char(value));                 
          end          
        else
          tm = char(34) + value(1) + char(34);
          for k = 2:length(value)
            tm = tm +  ", " + char(34) + value(k) + char(34);            
          end
          if length(value)>1
            tm = "[" + tm + "]";
          end
          fprintf(fid, '%s = %s;\n', key, char(tm));                
        end
      elseif iscell(value)
        if isempty(value)
          fprintf(fid, '%s = [];\n', key);       
        else
          fprintf(fid, '%s = <cell, size: [%s]>;\n', key, num2str(size(value)));
        end
      elseif isstruct(value)
          fprintf(fid, '%s = <struct with %d field(s)>;\n', key, numel(fieldnames(value)));
      else
          fprintf(fid, '%s = <unsupported type: %s>;\n', key, class(value));
      end
      
      if key == "platform" || key == "modelnumber" || key == "nstage" || key == "preconditioner" 
        fprintf(fid,'\n');
      end
      if key == "externalparam" || key == "nfb" 
        fprintf(fid,'\n');
      end
      
    end
end

end

  
  