/*
    writebinaryfiles.cpp

    This file provides functions for writing mesh and solution data to binary files, 
    as well as mesh partitioning utilities (with METIS support). It is designed for 
    use in high-performance scientific computing applications involving PDE solvers.

    Functions:

    - writesol: 
        Writes the solution data for a given mesh partition to a binary file. 
        Includes mesh coordinates, solution variables, and metadata.

    - writeDouble:
        Helper function to write a single double value to a binary file.

    - writeVectorAsDoubles:
        Helper function to write a vector of integers as doubles to a binary file.

    - max_of_vector_or_zero:
        Returns the maximum value in a vector of integers, or zero if the vector is empty.

    - writemesh:
        Writes mesh connectivity and partitioning information to a binary file. 
        Includes support for hybrid meshes and parallel partitioning.

    - partitionMesh (METIS only):
        Partitions the mesh into subdomains using METIS, filling element and node partition arrays.

    - writeBinaryFiles:
        Main entry point for writing all binary files required for a PDE simulation. 
        Handles directory creation, mesh building, partitioning, and file output for both 
        serial and parallel cases.

    Notes:
    - Uses C++ STL containers and file streams for binary I/O.
    - Requires METIS for mesh partitioning in parallel runs.
    - Relies on external types: PDE, Mesh, Master, DMD, Conn.
    - Assumes existence of utility functions: error, ensure_dir, make_path, buildMesh, buildConn, 
      build_dmdhdg, build_dmdldg, initializeDMD, freeCharArray, select_columns, mke2e, writepde, writemaster.

    Compilation:
    - Define HAVE_METIS to enable METIS-based partitioning.
    - Include this header in the main application to enable binary file output for mesh and solution data.
*/

#ifndef __WRITEBINARYFILES
#define __WRITEBINARYFILES

void writesol(const PDE& pde, const Mesh& mesh, const Master& master, const Conn& conn, const vector<int>& elempart, const vector<int>& facepartpts, const std::string& filename) 
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) error("Error opening file: " + filename);
    
    std::vector<double> ndims(12, 0);
    ndims[0] = elempart.size();
    ndims[1] = std::accumulate(facepartpts.begin(), facepartpts.end(), 0);
    ndims[2] = mesh.nfe;
    ndims[3] = master.npe;
    ndims[4] = master.npf;
    ndims[5] = pde.nc;
    ndims[6] = pde.ncu;
    ndims[7] = pde.ncq;
    ndims[8] = pde.ncw;
    ndims[9] = pde.ncv;
    ndims[10] = pde.nch;
    ndims[11] = pde.ncx;

    int ncudg = 0;

    std::vector<double> nsize(20, 0);
    nsize[0] = static_cast<double>(ndims.size());
    if (!mesh.xdg.empty()) {
        nsize[1] = static_cast<double>(master.npe*mesh.dim*elempart.size());        
    }
    if (!mesh.udg.empty()) {
        ncudg = mesh.udg.size()/(mesh.ne*master.npe);
        if (ncudg == pde.ncu || ncudg == pde.nc) 
            nsize[2] = static_cast<double>(master.npe*ncudg*elempart.size());
        else
            error("size of udg is incorrect.");
    }
    if (!mesh.vdg.empty()) {
        nsize[3] = static_cast<double>(master.npe*pde.ncv*elempart.size());
    }
    if (!mesh.wdg.empty()) {
        nsize[4] = static_cast<double>(master.npe*pde.ncw*elempart.size());
    }
    if (!mesh.uhat.empty()) {
        nsize[5] = static_cast<double>(mesh.uhat.size());
    }
    if (!conn.cgnodes.empty()) {
        nsize[6] = static_cast<double>(conn.cgnodes.size());
    }

    // Helper to write vector as binary double
    auto writeDoubleVector = [&](const std::vector<double>& vec) {
        file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));
    };

    double len = static_cast<double>(nsize.size());
    file.write(reinterpret_cast<const char*>(&len), sizeof(double));
    writeDoubleVector(nsize);
    writeDoubleVector(ndims);

    vector<double> tm;
    if (!mesh.xdg.empty()) {
        tm.resize(master.npe*mesh.dim*elempart.size());
        select_columns(tm.data(), mesh.xdg.data(), elempart.data(), master.npe*mesh.dim, elempart.size());
        writeDoubleVector(tm);
    }
    if (!mesh.udg.empty()) {
        tm.resize(master.npe*ncudg*elempart.size());
        select_columns(tm.data(), mesh.udg.data(), elempart.data(), master.npe*ncudg, elempart.size());
        writeDoubleVector(tm);    
    }
    if (!mesh.vdg.empty()) {
        tm.resize(master.npe*pde.ncv*elempart.size());
        select_columns(tm.data(), mesh.vdg.data(), elempart.data(), master.npe*pde.ncv, elempart.size());
        writeDoubleVector(tm);                    
    }
    if (!mesh.wdg.empty()) {
        tm.resize(master.npe*pde.ncw*elempart.size());
        select_columns(tm.data(), mesh.wdg.data(), elempart.data(), master.npe*pde.ncw, elempart.size());
        writeDoubleVector(tm);                    
    }
    if (!mesh.uhat.empty()) {
        writeDoubleVector(mesh.uhat);
    }
    if (!conn.cgnodes.empty()) {
        writeDoubleVector(conn.cgnodes);
    }
    
    file.close();
    
    std::cout << "Finished writing initial solution to " + filename << std::endl;
}

inline void writeDouble(std::ofstream& out, double v) {
    out.write(reinterpret_cast<const char*>(&v), sizeof(double));
}

void writeVectorAsDoubles(std::ofstream& out, const vector<int>& a) {
    for (size_t i = 0; i < a.size(); ++i) {
        double v = static_cast<double>(a[i]);
        out.write(reinterpret_cast<const char*>(&v), sizeof(double));
    }
}

inline int max_of_vector_or_zero(const std::vector<int>& v) {
    if (v.empty()) return 0;
    return *std::max_element(v.begin(), v.end());
}

void writemesh(const PDE& pde,
               const Mesh& mesh,
               const Master& master,               
               const DMD& dmd,
               const Conn& conn, 
               const std::string& filename)
{    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) error("Error opening file: " + filename);
    
    std::vector<int> ndims(20, 0);
    ndims[0] = (mesh.dim);
    ndims[1] = (dmd.elempart.size());
    ndims[2] = std::accumulate(conn.facepartpts.begin(), conn.facepartpts.end(), 0);
    ndims[3] = (max_of_vector_or_zero(mesh.t)+1);
    ndims[4] = (mesh.nfe);
    ndims[5] = (conn.nbe);
    ndims[6] = (conn.neb);
    ndims[7] = (conn.nbf);
    ndims[8] = (conn.nfb);
    
    std::vector<int> nsize(50, 0);
    nsize[0]  = 20; 
    nsize[1]  = (conn.facecon.size());
    nsize[2]  = (conn.eblks.size());
    nsize[3]  = (conn.fblks.size());
    if (pde.mpiprocs>1) {
      nsize[4]  = (dmd.nbsd.size());
      nsize[5]  = (dmd.elemsend.size());
      nsize[6]  = (dmd.elemrecv.size());
      nsize[7]  = (dmd.elemsendpts.size());
      nsize[8]  = (dmd.elemrecvpts.size());
    }
    nsize[9]  = (dmd.elempart.size());
    nsize[10] = (dmd.elempartpts.size());    
    nsize[11] = (conn.cgelcon.size());
    nsize[12] = (conn.rowent2elem.size());
    nsize[13] = (conn.cgent2dgent.size());
    nsize[14] = (conn.colent2elem.size());
    nsize[15] = (conn.rowe2f1.size());
    nsize[16] = (conn.cole2f1.size());
    nsize[17] = (conn.ent2ind1.size());
    nsize[18] = (conn.rowe2f2.size());
    nsize[19] = (conn.cole2f2.size());
    nsize[20] = (conn.ent2ind2.size());    
    nsize[21] = (conn.f2t.size());
    nsize[22] = (conn.elemcon.size());
    
    nsize[23] = (master.perm.size());
    nsize[24] = (conn.bf.size());
    nsize[25] = (mesh.cartGridPart.size());
    
    int ne = dmd.elempart.size();
    int nve = mesh.nve;
    nsize[26] = ne * nve;
    nsize[27] = mesh.boundaryConditions.size();
    nsize[28] = dmd.intepartpts.size();

    writeDouble(out, static_cast<double>(nsize.size())); 
    writeVectorAsDoubles(out, nsize);
    writeVectorAsDoubles(out, ndims);

    writeVectorAsDoubles(out, conn.facecon); // 
    writeVectorAsDoubles(out, conn.eblks);   //
    writeVectorAsDoubles(out, conn.fblks);   //  

    if (pde.mpiprocs>1) {      
      int n1 = dmd.elemsend.size();
      int n2 = dmd.elemrecv.size();
      vector<int> elemsend(n1); 
      vector<int> elemrecv(n2);
      for (int i=0; i<n1; i++) elemsend[i] = dmd.elemsend[i][1];
      for (int i=0; i<n2; i++) elemrecv[i] = dmd.elemrecv[i][1];
      writeVectorAsDoubles(out, dmd.nbsd);
      writeVectorAsDoubles(out, elemsend);
      writeVectorAsDoubles(out, elemrecv);
      writeVectorAsDoubles(out, dmd.elemsendpts);
      writeVectorAsDoubles(out, dmd.elemrecvpts);
    }
    writeVectorAsDoubles(out, dmd.elempart);
    writeVectorAsDoubles(out, dmd.elempartpts);

    writeVectorAsDoubles(out, conn.cgelcon);     //
    writeVectorAsDoubles(out, conn.rowent2elem); //
    writeVectorAsDoubles(out, conn.cgent2dgent); //
    writeVectorAsDoubles(out, conn.colent2elem); //
    writeVectorAsDoubles(out, conn.rowe2f1);     //
    writeVectorAsDoubles(out, conn.cole2f1);     //
    writeVectorAsDoubles(out, conn.ent2ind1);    //
    writeVectorAsDoubles(out, conn.rowe2f2);     //
    writeVectorAsDoubles(out, conn.cole2f2);     //
    writeVectorAsDoubles(out, conn.ent2ind2);    //

    if (pde.hybrid > 0) {
        writeVectorAsDoubles(out, conn.f2t);     //
        writeVectorAsDoubles(out, conn.elemcon); //
        writeVectorAsDoubles(out, master.perm);  
        writeVectorAsDoubles(out, conn.bf);      
        writeVectorAsDoubles(out, mesh.cartGridPart);
    }

    vector<int> ti(nve*ne); 
    select_columns(ti.data(), mesh.t.data(), dmd.elempart.data(), nve, ne); 
    writeVectorAsDoubles(out, ti);
    writeVectorAsDoubles(out, mesh.boundaryConditions);
    writeVectorAsDoubles(out, dmd.intepartpts);

    out.close();    
    std::cout << "Finished writing mesh to " + filename << std::endl;
}

#ifdef HAVE_METIS
void partitionMesh(std::vector<int>& epart, std::vector<int>& npart, std::vector<int>& eind, 
         int ne, int np, int nve, int nvf, int nparts) 
{ 
    std::vector<int> eptr(ne+1);
    eptr[0] = 0;
    for (int i=0; i<ne; i++) eptr[i+1] = eptr[i] + nve;
      
    // Resize output vectors
    epart.resize(ne);
    npart.resize(np);
  
    // METIS options
    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;  
    
    int objval = 0;
    int status = METIS_PartMeshDual(&ne, &np, eptr.data(), eind.data(),
                                    NULL, NULL, &nvf, &nparts,
                                    NULL, options, &objval,
                                    epart.data(), npart.data());
    
    if (status != METIS_OK) {
        std::cerr << "METIS_PartMeshDual failed with status " << status << ".\n";
    }  else {
        std::cout << "Finished partitioning mesh using METIS" << std::endl;
    }    
}
#endif

void writeBinaryFiles(PDE& pde, Mesh& mesh, const Master& master, const ParsedSpec& spec) 
{
    bool callbuildConn = false;

    ensure_dir(pde.datainpath);
    ensure_dir(pde.dataoutpath);
    
    for (const auto& vec : spec.vectors) {
        const std::string& name = vec.first;
        int size = vec.second;
        //std::cout<<name<<" : "<<size<<std::endl;
        if (name == "uhat") pde.ncu = size;
        if (name == "v") pde.ncv = size;
        if (name == "w") pde.ncw = size;
        if (name == "uq") pde.nc = size;        
    }

    for (int i=0; i<spec.functions.size(); i++) {
        //std::cout<<spec.functions[i].name<<" : "<<spec.functions[i].outputsize<<std::endl;
        if (spec.functions[i].name == "VisScalars") pde.nsca = spec.functions[i].outputsize;
        if (spec.functions[i].name == "VisVectors") pde.nvec = spec.functions[i].outputsize/pde.nd;
        if (spec.functions[i].name == "VisTensors") pde.nten = spec.functions[i].outputsize/(pde.nd*pde.nd);
        if (spec.functions[i].name == "QoIboundary") pde.nsurf = spec.functions[i].outputsize;
        if (spec.functions[i].name == "QoIvolume") pde.nvqoi = spec.functions[i].outputsize;
    }
    
    writepde(pde, make_path(pde.datainpath, "app.bin"));
    writemaster(master, make_path(pde.datainpath, "master.bin"));    

    if (pde.writemeshsol == 1) {
        buildMesh(mesh, pde, master);
            
        if (pde.mpiprocs>1) {
        
            if ((pde.partitionfile == "") || (mesh.elem2cpu.size() == 0)) {
#ifdef HAVE_METIS          
                vector<int> node2cpu;
                partitionMesh(mesh.elem2cpu, node2cpu, mesh.t, mesh.ne, mesh.np, mesh.nve, mesh.nvf, pde.mpiprocs);
                node2cpu.resize(0);
                for (int i=0; i<mesh.ne; i++) mesh.elem2cpu[i] += 1;        
#else
                error("mpiprocs > 1 requires a mesh partition array. \nPlease include the required mesh partition array in a binary file\nand set partitionfile to the name of the file.");      
#endif                  
            }

            for (int i=0; i<mesh.ne; i++) mesh.elem2cpu[i] -= 1;                    
            mesh.t2t.resize(mesh.nfe*mesh.ne);
            mesh.nf = mke2e_hash(mesh.t2t.data(), mesh.t.data(), mesh.localfaces.data(), mesh.ne, mesh.nve, mesh.nvf, mesh.nfe);
            
        //       if (pde.debugmode==1) {
        //         writearray2file(pde.datapath + "/t2t.bin", mesh.t2t.data(), mesh.t2t.size());
        //         writearray2file(pde.datapath + "/elem2cpu.bin", mesh.elem2cpu.data(), mesh.elem2cpu.size());
        //       }
            
            //cout<<mesh.nf<<", "<<pde.coupledinterface<<": here 1"<<endl;
            vector<DMD> dmd(pde.mpiprocs);            
            if (pde.hybrid==1)         
                build_dmdhdg(dmd, mesh.t2t.data(), mesh.elem2cpu.data(), mesh.inte.data(), mesh.nfe, mesh.ne, pde);
            else build_dmdldg(dmd, mesh.t2t.data(), mesh.elem2cpu.data(), mesh.nfe, mesh.ne, pde);                
            
            for (int n=0; n<pde.mpiprocs; n++) {
                Conn conn;          
                if (callbuildConn)  buildConn(conn, pde, mesh, master, dmd[n]);
                else {
                    int ne = dmd[n].elempart.size();
                    conn.bf.resize(mesh.nfe * ne); 
                    int* fi = (int*)malloc(mesh.nfe * ne * sizeof(int));      
                    select_columns(fi, mesh.f.data(), dmd[n].elempart.data(), mesh.nfe, ne);       
                    apply_bcm(conn.bf.data(), fi, mesh.boundaryConditions.data(), mesh.nfe*ne, mesh.nbcm);                              
                    CPUFREE(fi);
                }
                
                writesol(pde, mesh, master, conn, dmd[n].elempart, conn.facepartpts, make_path(pde.datainpath, "sol" + std::to_string(n+1) + ".bin"));    
                writemesh(pde, mesh, master, dmd[n], conn, make_path(pde.datainpath, "mesh" + std::to_string(n+1) + ".bin"));
            }               
        } else {
            DMD dmd = initializeDMD(pde, mesh);    
            Conn conn;                
            if (callbuildConn) buildConn(conn, pde, mesh, master, dmd);
            else {
                int ne = dmd.elempart.size();
                conn.bf.resize(mesh.nfe * ne); 
                int* fi = (int*)malloc(mesh.nfe * ne * sizeof(int));      
                select_columns(fi, mesh.f.data(), dmd.elempart.data(), mesh.nfe, ne);       
                apply_bcm(conn.bf.data(), fi, mesh.boundaryConditions.data(), mesh.nfe*ne, mesh.nbcm);                              
                CPUFREE(fi);
            }
    
            writesol(pde, mesh, master, conn, dmd.elempart, conn.facepartpts, make_path(pde.datainpath , "sol.bin"));    
            writemesh(pde, mesh, master, dmd, conn, make_path(pde.datainpath, "mesh.bin"));      
        }    
    }

    freeCharArray(mesh.boundaryExprs, mesh.nbndexpr);
    freeCharArray(mesh.curvedBoundaryExprs, mesh.nbndexpr);
    freeCharArray(mesh.periodicExprs1, mesh.nprdexpr*mesh.nprdcom);
    freeCharArray(mesh.periodicExprs2, mesh.nprdexpr*mesh.nprdcom);                
    
    std::cout << "Finished writeBinaryFiles.\n";
}

#endif
