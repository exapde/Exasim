# Helpful cases: in the examples directory:
- Poisson/Poisson2D
- NavierStokes/naca2d
- NavierStokes/hypersoniccylinder_mach8
  - 2D steady state high-speed flow
  - Read in inputs, run full continuation solver
- NavierStokes/nshtmach8
  - Probably most relevant

# General notes: 
- pdeapp.m is usually where the code is run
- fluxes, etc. defined in pdemodel, usuaully
- `exasim` function generates code, compiles, runs. Can be more convenient to break up this API and call it in separate parts. For example
    
    ```
        pde = setcompilers(pde);       
        
        % generate input files and store them in datain folder
        [pde,mesh,master,dmd] = preprocessing(pde,mesh);

        % generate source codes and store them in app folder
        pde.gencode = 1 % or 0
        if pde.gencode==1
        kkgencode(pde);
        end
        
        compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
        
        runstr = runcode(pde, 1); % run C++ code

        % get solution from output files in dataout folder
        sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
    ```

# Notes on coupling
- Currently, the way we transfer solution fields into an Exasim solver is with the mesh.vdg field. 
- Right now, this is used in MATLAB. Could be read from a file. 
- Easiest way for now: read in data into Matlab, put in mesh.vdg, etc.
- Most relevant file may be `solutiontransfer_ht2ns.m`
  - What does it do? Transfers wall temperature from solid to a fluid solver
- Also useful: `solutiontransfer_ns2ht.m`
  - What does it do? Transfer heat flux from fluid to solid tomain

# TODOs
- How should Summit solution info be saved for Exasim to read it?  
  - First workflow will be to save Summit output to a file and read to Exasim. Summit needs to know how to save it 
  - Summit saves vtk...can we read that directly? 