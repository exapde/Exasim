function saveStructToTextFile(s, filename)
% saveStructToTextFile(s, filename)
% Writes struct fields and values into a text file.

fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file %s for writing.', filename);
end

%fprintf(fid, 'Struct saved on: %s\n\n', datestr(now));
fields = fieldnames(s);
requiredKeys = ["exasimpath", "datapath", "model", "modelfile", "meshfile", "discretization",...
        "platform", "mpiprocs", "porder", "pgauss", "physicsparam",...
        "tau", "boundaryconditions", "boundaryexpressions"];

for i = 1:length(fields)
    key = fields{i};
    value = s.(key);
    
    match = 0;
    for j = 1:length(requiredKeys)
      if (key == requiredKeys(j)) 
        match = 1;
      end
    end
      
    if match == 1
      if isnumeric(value) || islogical(value)
          fprintf(fid, '%s = %s;\n', key, mat2str(value));
      elseif ischar(value)
          fprintf(fid, '%s = ''%s'';\n', key, value);
      elseif isstring(value)
          fprintf(fid, '%s = \"%s\";\n', key, char(value));
      elseif iscell(value)
          fprintf(fid, '%s = <cell, size: [%s]>;\n', key, num2str(size(value)));
      elseif isstruct(value)
          fprintf(fid, '%s = <struct with %d field(s)>;\n', key, numel(fieldnames(value)));
      else
          fprintf(fid, '%s = <unsupported type: %s>;\n', key, class(value));
      end
    end
end

fclose(fid);
fprintf('Struct saved to %s\n', filename);
end
