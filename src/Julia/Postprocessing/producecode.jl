function  producecode(app,mesh)

# search compilers and set options
app = Gencode.setcompilers(app);

# generate input files and store them in datain folder
app, mesh, master, dmd = preprocessing(app,mesh);

# generate source codes and store them in app folder
Gencode.gencode(app);

# compile source codes to build an executable file and store it in app folder
compilerstr = Gencode.compilecode(app);

return compilerstr,app,mesh,master,dmd

end
