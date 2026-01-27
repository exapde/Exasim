function addBuiltinModel(app)

if app.builtinmodelID > 0
  disp("This PDE model is added to the builtin model library ...");
  
  id = app.builtinmodelID;
  cdir = pwd();
  fdir = app.exasimpath + "/frontends/Matlab/Gencode/";
  mdir = app.backendpath + "/Model/BuiltIn/";  
  ndir = mdir + num2str(id);

  text = fileread(char(cdir + "/" + app.modelfile + ".m"));    
  fid = fopen(mdir + "pdemodel" + num2str(id) + ".m", 'w');  
  t = string(text);          % ensure scalar string
  t = strrep(t, "pdemodel", "pdemodel" + num2str(id));
  t = replace(t, newline, sprintf('\n'));  % normalize newlines
  fwrite(fid, char(t), 'char');
  fclose(fid);

  text = fileread(char(fdir + "genmodel.m"));    
  fid = fopen(mdir + "genmodel" + num2str(id) + ".m", 'w');  
  t = string(text);          % ensure scalar string
  t = strrep(t, "pde.builtinmodelID = 0", "pde.builtinmodelID = " + num2str(id));
  t = strrep(t, "pde.nd = 0", "pde.nd = " + num2str(app.nd));
  t = strrep(t, "pde.ncu = 0", "pde.ncu = " + num2str(app.ncu));
  t = strrep(t, "pde.ncq = 0", "pde.ncq = " + num2str(app.ncq));
  t = strrep(t, "pde.ncw = 0", "pde.ncw = " + num2str(app.ncw));
  t = strrep(t, "pde.ncv = 0", "pde.ncv = " + num2str(app.nco));
  t = strrep(t, "pde.ntau = 0", "pde.ntau = " + num2str(numel(app.tau)));
  t = strrep(t, "pde.nmu = 0", "pde.nmu = " + num2str(numel(app.physicsparam)));
  t = strrep(t, "pde.neta = 0", "pde.neta = " + num2str(numel(app.externalparam)));  
  t = replace(t, newline, sprintf('\n'));  % normalize newlines
  fwrite(fid, char(t), 'char');
  fclose(fid);
  
  run(mdir + "genmodel" + num2str(id) + ".m");  
end

