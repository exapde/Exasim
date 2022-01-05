import Preprocessing, Gencode

def producecode(pde,mesh):

    # search compilers and set options
    pde = Gencode.setcompilers(pde);

    # generate input files and store them in datain folder
    pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

    # generate source codes and store them in app folder
    Gencode.gencode(pde);

    # compile source codes to build an executable file and store it in app folder
    compilerstr = Gencode.compilecode(pde);

    return compilerstr,pde,mesh,master,dmd
