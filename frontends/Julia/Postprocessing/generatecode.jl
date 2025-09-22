function  generatecode(app,mesh)

    # search compilers and set options
    app = Gencode.setcompilers(app);

    # generate input files and store them in datain folder
    app, mesh, master, dmd = preprocessing(app,mesh);
    
    # generate source codes and store them in app folder
    Gencode.gencode(app);
        
    return app,mesh,master,dmd
    
end
    