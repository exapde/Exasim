function [match, mesh1, mesh2] = checkmeshbin(file1, file2)

a = readbin(file1);
b = readbin(file2);
match = 0;
if numel(a) == numel(b)
  fprintf("check mesh: %g\n",max(abs(a-b)));
  if max(abs(a-b)) <= 1e-10 
    match = 1;
  end
end

mesh1 = readmeshbin(file1);
mesh2 = readmeshbin(file2);

fields = fieldnames(mesh1);
for i = 1:length(fields)
    field_name = fields{i};
    a = mesh1.(field_name); 
    b = mesh2.(field_name); 
    a = a(:); b = b(:); s = field_name;
    if numel(a) == numel(b)      
      fprintf("Check " + s + ":  (%g, %g)\n", max(abs(a(:)-b(:))), max(abs(1+a(:)-b(:))));      
    else
      fprintf(s + ": sizes do not match (%d, %d)\n", numel(a), numel(b));
    end
end
