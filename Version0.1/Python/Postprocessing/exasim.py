import Preprocessing, Gencode
from fetchsolution import fetchsolution

def exasim(pde,mesh):

    # search compilers and set options
    pde = Gencode.setcompilers(pde);

    # generate input files and store them in datain folder
    pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

    # generate source codes and store them in app folder
    Gencode.gencode(pde);

    # compile source codes to build an executable file and store it in app folder
    compilerstr = Gencode.compilecode(pde);

    # run executable file to compute solution and store it in dataout folder
    runstr = Gencode.runcode(pde);

    # get solution from output files in dataout folder
    pde['vistime'] = [];
    sol = fetchsolution(pde,master,dmd);

    return sol,pde,mesh,master,dmd,compilerstr,runstr
