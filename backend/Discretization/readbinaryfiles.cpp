/*
    readbinaryfiles.cpp

    This file provides functions for reading and writing binary files containing the data structures
    used in the Exasim backend for high-order finite element simulations. The main structures handled
    are appstruct, masterstruct, meshstruct, and solstruct, which store application parameters, master
    element data, mesh connectivity, and solution data, respectively.

    Functions:

    - void readappstruct(string filename, appstruct &app)
        Reads application parameters from a binary file into an appstruct object. Handles initialization
        of physics and solver parameters, and optionally initializes Mutation++ mixture options if enabled.

    - void writeappstruct(string filename, appstruct &app)
        Writes the contents of an appstruct object to a binary file.

    - void readmasterstruct(string filename, masterstruct &master)
        Reads master element data from a binary file into a masterstruct object.

    - void writemasterstruct(string filename, masterstruct &master)
        Writes the contents of a masterstruct object to a binary file.

    - void readmeshstruct(string filename, meshstruct &mesh)
        Reads mesh connectivity and partitioning data from a binary file into a meshstruct object.

    - void writemeshstruct(string filename, meshstruct &mesh)
        Writes the contents of a meshstruct object to a binary file.

    - void readsolstruct(string filename, solstruct &sol)
        Reads solution data (coordinates, unknowns, auxiliary variables) from a binary file into a solstruct object.

    - void writesolstruct(string filename, solstruct &sol)
        Writes the contents of a solstruct object to a binary file.

    - void readsolstruct(string filename, solstruct &sol, appstruct &app, masterstruct &master, meshstruct &mesh, Int mpirank)
        Reads solution data from a binary file, with additional logic for initializing solution arrays
        based on the application, master, and mesh structures, and MPI rank.

    - void readInput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string filein, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank)
        Reads all input data structures from binary files, handling both serial and parallel (MPI) cases.

    - void writeOutput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string fileout, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank)
        Writes all output data structures to binary files, handling both serial and parallel (MPI) cases.

    Dependencies:
        - Standard C++ libraries: <string>, <fstream>
        - Optional: Mutation++ library (if HAVE_MPP is defined)
        - Custom types: appstruct, masterstruct, meshstruct, solstruct, Int, dstype
        - Utility functions: readiarrayfromdouble, readarray, writeiarraytodouble, writearray, error, NumberToString, cpuArraySetValue, cpuArrayInsert, cpuInituDriver, cpuInitudgDriver, cpuInitodgDriver, cpuInitwdgDriver, CPUFREE

    Notes:
        - All file operations are performed in binary mode for efficiency and compatibility.
        - Error handling is performed via the error() function.
        - MPI and OpenMP support is handled via function arguments and conditional logic.
        - Some functions contain conditional compilation for Mutation++ and Enzyme support.
*/
#ifndef __READBINARYFILES
#define __READBINARYFILES

#ifdef HAVE_MPP
#include <mutation++.h>
#endif

void readappstruct(string filename, appstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) 
        error("Unable to open file " + filename);
       
    //printf("Read app struct from files...\n");   
    
    /* Read data to app structure */            
    app.lsize = readiarrayfromdouble(in, 1);
    app.nsize = readiarrayfromdouble(in, app.lsize[0]);
    app.ndims = readiarrayfromdouble(in, app.nsize[0]);
    app.flag = readiarrayfromdouble(in, app.nsize[1]);
    app.problem = readiarrayfromdouble(in, app.nsize[2]);
    readarray(in, &app.uinf, app.nsize[3]);
    readarray(in, &app.dt, app.nsize[4]);                
    readarray(in, &app.factor, app.nsize[5]);       
    readarray(in, &app.physicsparam, app.nsize[6]);       
    readarray(in, &app.solversparam, app.nsize[7]);   
    readarray(in, &app.tau, app.nsize[8]);   
    readarray(in, &app.stgdata, app.nsize[9]);   
    readarray(in, &app.stgparam, app.nsize[10]);   
    app.stgib = readiarrayfromdouble(in, app.nsize[11]);
    app.vindx = readiarrayfromdouble(in, app.nsize[12]);
    readarray(in, &app.dae_dt, app.nsize[13]);   
    app.interfacefluxmap = readiarrayfromdouble(in, app.nsize[14]);
    readarray(in, &app.avparam, app.nsize[15]);   
    
    app.szflag = app.nsize[1];
    app.szproblem = app.nsize[2];
    app.szuinf = app.nsize[3];
    app.szdt = app.nsize[4];
    app.szfactor = app.nsize[5];
    app.szphysicsparam = app.nsize[6];
    app.szsolversparam = app.nsize[7];
    app.sztau = app.nsize[8];
    app.szstgdata = app.nsize[9];
    app.szstgparam = app.nsize[10];
    app.szstgib = app.nsize[11];
    app.szvindx = app.nsize[12];
    app.szdae_dt = app.nsize[13];
    app.szinterfacefluxmap = app.nsize[14];
    app.szavparam = app.nsize[15];

    #ifdef HAVE_MPP
        char a[50];
        in.getline(a, 50, 'X');
        string MixtureName = a;
        
        char b[50];
        in.getline(b, 50, 'X');
        string StateModel = b;

        char c[50];
        in.getline(c, 50, 'X');
        string ThermoDB = c;

        Mutation::MixtureOptions opts(MixtureName);
        opts.setStateModel(StateModel);
        opts.setThermodynamicDatabase(ThermoDB);
        opts.setViscosityAlgorithm("Gupta-Yos");
        opts.setThermalConductivityAlgorithm("Wilke");
        
        app.mix = new Mutation::Mixture(opts);
        printf("Mutation mixture initialized\n");
    #endif 

    Int i, ncu, ncq, ncw;
    ncu = app.ndims[6];// number of compoments of (u)
    ncq = app.ndims[7];// number of compoments of (q)
    ncw = app.ndims[13];// number of compoments of (w)
    
    if (ncu>0) {
        app.fc_u = (dstype*) malloc(sizeof(dstype)*ncu);
        app.dtcoef_u = (dstype*) malloc(sizeof(dstype)*ncu);
        for (i=0; i<ncu; i++) {
            app.fc_u[i] = 1.0;     //app.factor[i];
            app.dtcoef_u[i] = 1.0; //app.factor[i];
        }
        app.szfc_u = ncu;
        app.szdtcoef_u = ncu;
    }        
    if (ncq>0) {
        app.fc_q = (dstype*) malloc(sizeof(dstype)*ncq);
        app.dtcoef_q = (dstype*) malloc(sizeof(dstype)*ncq);
        for (i=0; i<ncq; i++) {
            app.fc_q[i] = 1.0;     //app.factor[ncu+i];
            app.dtcoef_q[i] = 1.0; //app.factor[ncu+i];
        }
        app.szfc_q = ncu;
        app.szdtcoef_q = ncu;
    }        
    if (ncw>0) {
        app.fc_w = (dstype*) malloc(sizeof(dstype)*ncw);
        app.dtcoef_w = (dstype*) malloc(sizeof(dstype)*ncw);
        for (i=0; i<ncw; i++) {
            app.fc_w[i] = 1.0;    // app.factor[ncu+ncq+i];
            app.dtcoef_w[i] =1.0; // app.factor[ncu+ncq+i];
        }
        app.szfc_w = ncw;
        app.szdtcoef_w = ncw;
    }                
    //app.time = app.factor[ncu+ncq+ncp];

    // Close file:
    in.close();
}

void writeappstruct(string filename, appstruct &app)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);

    if (!out) {
        error("Unable to open file " + filename);
    }
       
    printf("Write app struct into file...\n");  
    
    /* write app structure to data file */      
    writeiarraytodouble(out, app.lsize, 1);
    writeiarraytodouble(out, app.nsize, app.lsize[0]);
    writeiarraytodouble(out, app.ndims, app.nsize[0]);
    writeiarraytodouble(out, app.flag, app.nsize[1]);
    writeiarraytodouble(out, app.problem, app.nsize[2]);
    writearray(out, app.uinf, app.nsize[3]);
    writearray(out, app.dt, app.nsize[4]);        
    writearray(out, app.factor, app.nsize[5]);       
    writearray(out, app.physicsparam, app.nsize[6]);       
    writearray(out, app.solversparam, app.nsize[7]);    
    writearray(out, app.tau, app.nsize[8]);      
    writearray(out, app.stgdata, app.nsize[9]);      
    writearray(out, app.stgparam, app.nsize[10]);      
    writeiarraytodouble(out, app.stgib, app.nsize[11]);
    writeiarraytodouble(out, app.vindx, app.nsize[12]);
    writearray(out, app.dae_dt, app.nsize[13]);        
    
    // Close file:
    out.close();
}

void readmasterstruct(string filename, masterstruct &master)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    //printf("Read master struct from files...\n");  
    
    /* Read data to app structure */        
    master.lsize = readiarrayfromdouble(in, 1);
    master.nsize = readiarrayfromdouble(in, master.lsize[0]);        
    master.ndims = readiarrayfromdouble(in, master.nsize[0]);
    readarray(in, &master.shapegt, master.nsize[1]);        
    readarray(in, &master.shapegw, master.nsize[2]);
    readarray(in, &master.shapfgt, master.nsize[3]);        
    readarray(in, &master.shapfgw, master.nsize[4]);    
    readarray(in, &master.shapent, master.nsize[5]);        
    readarray(in, &master.shapen, master.nsize[6]);
    readarray(in, &master.shapfnt, master.nsize[7]);        
    readarray(in, &master.shapfn, master.nsize[8]);
    readarray(in, &master.xpe, master.nsize[9]);        
    readarray(in, &master.gpe, master.nsize[10]);    
    readarray(in, &master.gwe, master.nsize[11]);        
    readarray(in, &master.xpf, master.nsize[12]);
    readarray(in, &master.gpf, master.nsize[13]);        
    readarray(in, &master.gwf, master.nsize[14]);                        
    readarray(in, &master.shap1dgt, master.nsize[15]);
    readarray(in, &master.shap1dgw, master.nsize[16]);
    readarray(in, &master.shap1dnt, master.nsize[17]);
    readarray(in, &master.shap1dnl, master.nsize[18]);
    readarray(in, &master.xp1d, master.nsize[19]);
    readarray(in, &master.gp1d, master.nsize[20]);
    readarray(in, &master.gw1d, master.nsize[21]);

    master.szshapegt = master.nsize[1];
    master.szshapegw = master.nsize[2];
    master.szshapfgt = master.nsize[3];
    master.szshapfgw = master.nsize[4];
    master.szshapent = master.nsize[5];
    master.szshapen = master.nsize[6];
    master.szshapfnt = master.nsize[7];
    master.szshapfn = master.nsize[8];
    master.szxpe = master.nsize[9];
    master.szgpe = master.nsize[10];
    master.szgwe = master.nsize[11];
    master.szxpf = master.nsize[12];
    master.szgpf = master.nsize[13];
    master.szgwf = master.nsize[14];
    master.szshap1dgt = master.nsize[15];
    master.szshap1dgw = master.nsize[16];
    master.szshap1dnt = master.nsize[17];
    master.szshap1dnl = master.nsize[18];
    master.szxp1d = master.nsize[19];
    master.szgp1d = master.nsize[20];
    master.szgw1d = master.nsize[21];

//     master.nd = master.ndims[0]; 
//     master.porder = master.ndims[3]; 
//     master.pgauss = master.ndims[4]; 
//     master.npe = master.ndims[5];
//     master.npf = master.ndims[6];
//     master.nge = master.ndims[7];
//     master.ngf = master.ndims[8];
    
    // Close file:
    in.close();
}


void writemasterstruct(string filename, masterstruct &master)
{
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);
    
    if (!out) {
        error("Unable to open file " + filename);
    }
        
    printf("Write master struct into files...\n");  
    
    writeiarraytodouble(out, master.lsize, 1);
    writeiarraytodouble(out, master.nsize, master.lsize[0]);
    writeiarraytodouble(out, master.ndims, master.nsize[0]);    
    writearray(out, master.shapegt, master.nsize[1]);        
    writearray(out, master.shapegw, master.nsize[2]);
    writearray(out, master.shapfgt, master.nsize[3]);        
    writearray(out, master.shapfgw, master.nsize[4]);    
    writearray(out, master.shapent, master.nsize[5]);        
    writearray(out, master.shapen, master.nsize[6]);
    writearray(out, master.shapfnt, master.nsize[7]);        
    writearray(out, master.shapfn, master.nsize[8]);
    writearray(out, master.xpe, master.nsize[9]);        
    writearray(out, master.gpe, master.nsize[10]);    
    writearray(out, master.gwe, master.nsize[11]);        
    writearray(out, master.xpf, master.nsize[12]);
    writearray(out, master.gpf, master.nsize[13]);        
    writearray(out, master.gwf, master.nsize[14]);    
    writearray(out, master.shap1dgt, master.nsize[15]);
    writearray(out, master.shap1dgw, master.nsize[16]);
    writearray(out, master.shap1dnt, master.nsize[17]);
    writearray(out, master.shap1dnl, master.nsize[18]);
    writearray(out, master.xp1d, master.nsize[19]);
    writearray(out, master.gp1d, master.nsize[20]);
    writearray(out, master.gw1d, master.nsize[21]);
    
    // Close file:
    out.close();
}


void readmeshstruct(string filename, meshstruct &mesh, solstruct &sol, appstruct &app, masterstruct &master, Int mpirank)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
        
    mesh.lsize = readiarrayfromdouble(in, 1);
    mesh.nsize = readiarrayfromdouble(in, mesh.lsize[0]);
    mesh.ndims = readiarrayfromdouble(in, mesh.nsize[0]);
    mesh.facecon = readiarrayfromdouble(in, mesh.nsize[1]); //  
    mesh.eblks = readiarrayfromdouble(in, mesh.nsize[2]);  //
    mesh.fblks = readiarrayfromdouble(in, mesh.nsize[3]);  //
    mesh.nbsd = readiarrayfromdouble(in, mesh.nsize[4]);
    mesh.elemsend = readiarrayfromdouble(in, mesh.nsize[5]);
    mesh.elemrecv = readiarrayfromdouble(in, mesh.nsize[6]);
    mesh.elemsendpts = readiarrayfromdouble(in, mesh.nsize[7]);
    mesh.elemrecvpts = readiarrayfromdouble(in, mesh.nsize[8]);
    mesh.elempart = readiarrayfromdouble(in, mesh.nsize[9]);
    mesh.elempartpts = readiarrayfromdouble(in, mesh.nsize[10]);    
    mesh.cgelcon = readiarrayfromdouble(in, mesh.nsize[11]); //
    mesh.rowent2elem = readiarrayfromdouble(in, mesh.nsize[12]); //
    mesh.cgent2dgent = readiarrayfromdouble(in, mesh.nsize[13]); //
    mesh.colent2elem = readiarrayfromdouble(in, mesh.nsize[14]); //
    mesh.rowe2f1 = readiarrayfromdouble(in, mesh.nsize[15]);     //
    mesh.cole2f1 = readiarrayfromdouble(in, mesh.nsize[16]);    //
    mesh.ent2ind1 = readiarrayfromdouble(in, mesh.nsize[17]);   //
    mesh.rowe2f2 = readiarrayfromdouble(in, mesh.nsize[18]);    //
    mesh.cole2f2 = readiarrayfromdouble(in, mesh.nsize[19]);    //
    mesh.ent2ind2 = readiarrayfromdouble(in, mesh.nsize[20]);   //
    mesh.f2e = readiarrayfromdouble(in, mesh.nsize[21]);        //
    mesh.elemcon = readiarrayfromdouble(in, mesh.nsize[22]);    //
    mesh.perm = readiarrayfromdouble(in, mesh.nsize[23]);
    mesh.bf = readiarrayfromdouble(in, mesh.nsize[24]);
    mesh.cartgridpart = readiarrayfromdouble(in, mesh.nsize[25]);
  
    mesh.szfacecon = mesh.nsize[1];
    mesh.szeblks = mesh.nsize[2];
    mesh.szfblks = mesh.nsize[3];
    mesh.sznbsd = mesh.nsize[4];
    mesh.szelemsend = mesh.nsize[5];
    mesh.szelemrecv = mesh.nsize[6];
    mesh.szelemsendpts = mesh.nsize[7];
    mesh.szelemrecvpts = mesh.nsize[8];
    mesh.szelempart = mesh.nsize[9];
    mesh.szelempartpts = mesh.nsize[10];
    mesh.szcgelcon = mesh.nsize[11];
    mesh.szrowent2elem = mesh.nsize[12];
    mesh.szcgent2dgent = mesh.nsize[13];
    mesh.szcolent2elem = mesh.nsize[14];
    mesh.szrowe2f1 = mesh.nsize[15];
    mesh.szcole2f1 = mesh.nsize[16];
    mesh.szent2ind1 = mesh.nsize[17];
    mesh.szrowe2f2 = mesh.nsize[18];
    mesh.szcole2f2 = mesh.nsize[19];
    mesh.szent2ind2 = mesh.nsize[20];
    mesh.szf2e = mesh.nsize[21];
    mesh.szelemcon = mesh.nsize[22];
    mesh.szperm = mesh.nsize[23];
    mesh.szbf = mesh.nsize[24];
    mesh.szcartgridpart = mesh.nsize[25];
        
    int *ti = readiarrayfromdouble(in, mesh.nsize[26]);
    mesh.boundaryConditions = readiarrayfromdouble(in, mesh.nsize[27]);
    mesh.intepartpts = readiarrayfromdouble(in, mesh.nsize[28]);

    //printf("%d %d %d\n", mesh.nsize[26], mesh.nsize[27], mesh.nsize[28]);
    //checkConn(mesh, sol, app, master, ti, boundaryConditions, intepartpts, mesh.nsize[27]);    
    if (mesh.nsize[26] > 0 && mesh.nsize[27] > 0 && mesh.szfacecon == 0) {
        if (mpirank==0) printf("Building element and face connectivities \n");         
        buildConn(mesh, sol, app, master, ti, mesh.boundaryConditions, mesh.intepartpts, mesh.nsize[27]);
        CPUFREE(ti);
    }
    
    mesh.faceperm = readiarrayfromdouble(in, mesh.nsize[39]);
    mesh.nbintf = readiarrayfromdouble(in, mesh.nsize[40]);
    mesh.facesend = readiarrayfromdouble(in, mesh.nsize[41]);
    mesh.facesendpts = readiarrayfromdouble(in, mesh.nsize[42]);
    mesh.facerecv = readiarrayfromdouble(in, mesh.nsize[43]);
    mesh.facerecvpts = readiarrayfromdouble(in, mesh.nsize[44]);    
    mesh.szfaceperm = mesh.nsize[39];
    mesh.sznbintf = mesh.nsize[40];
    mesh.szfacesend = mesh.nsize[41];
    mesh.szfacesendpts = mesh.nsize[42];
    mesh.szfacerecv = mesh.nsize[43];
    mesh.szfacerecvpts = mesh.nsize[44];
        
    // Close file:
    in.close();            
}


void writemeshstruct(string filename, meshstruct &mesh)
{    
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);
    
    if (!out) {
        error("Unable to open file " + filename);
    }

    printf("Write mesh struct into files...\n");  
    
    writeiarraytodouble(out, mesh.lsize, 1);
    writeiarraytodouble(out, mesh.nsize, mesh.lsize[0]);    
    writeiarraytodouble(out, mesh.ndims, mesh.nsize[0]);
    writeiarraytodouble(out, mesh.facecon, mesh.nsize[1]);
    writeiarraytodouble(out, mesh.eblks, mesh.nsize[2]);
    writeiarraytodouble(out, mesh.fblks, mesh.nsize[3]);
    writeiarraytodouble(out, mesh.nbsd, mesh.nsize[4]);
    writeiarraytodouble(out, mesh.elemsend, mesh.nsize[5]);
    writeiarraytodouble(out, mesh.elemrecv, mesh.nsize[6]);
    writeiarraytodouble(out, mesh.elemsendpts, mesh.nsize[7]);
    writeiarraytodouble(out, mesh.elemrecvpts, mesh.nsize[8]);
    writeiarraytodouble(out, mesh.elempart, mesh.nsize[9]);
    writeiarraytodouble(out, mesh.elempartpts, mesh.nsize[10]);
    writeiarraytodouble(out, mesh.cgelcon, mesh.nsize[11]);
    writeiarraytodouble(out, mesh.rowent2elem, mesh.nsize[12]);
    writeiarraytodouble(out, mesh.cgent2dgent, mesh.nsize[13]);
    writeiarraytodouble(out, mesh.colent2elem, mesh.nsize[14]);
    writeiarraytodouble(out, mesh.rowe2f1, mesh.nsize[15]);
    writeiarraytodouble(out, mesh.cole2f1, mesh.nsize[16]);
    writeiarraytodouble(out, mesh.ent2ind1, mesh.nsize[17]);
    writeiarraytodouble(out, mesh.rowe2f2, mesh.nsize[18]);
    writeiarraytodouble(out, mesh.cole2f2, mesh.nsize[19]);
    writeiarraytodouble(out, mesh.ent2ind2, mesh.nsize[20]);
    
    // Close file:
    out.close();    
}

void readsolstruct(string filename, solstruct &sol)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    //printf("Read sol struct from files...\n");  

    sol.lsize = readiarrayfromdouble(in, 1);
    sol.nsize = readiarrayfromdouble(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdouble(in, sol.nsize[0]);    
    readarray(in, &sol.xdg, sol.nsize[1]);
    readarray(in, &sol.udg, sol.nsize[2]);    
    readarray(in, &sol.odg, sol.nsize[3]);  
    readarray(in, &sol.wdg, sol.nsize[4]);  
    
    // Close file:
    in.close();            
}

void writesolstruct(string filename, solstruct &sol)
{    
    // Open file to write
    ofstream out(filename.c_str(), ios::out | ios::binary);
    
    if (!out) {
        error("Unable to open file " + filename);
    }

    printf("Write sol struct into files...\n");  
    
    writeiarraytodouble(out, sol.lsize, 1);
    writeiarraytodouble(out, sol.nsize, sol.lsize[0]);            
    writeiarraytodouble(out, sol.ndims, sol.nsize[0]);
    writearray(out, sol.xdg, sol.nsize[1]);
    writearray(out, sol.udg, sol.nsize[2]);    
    writearray(out, sol.odg, sol.nsize[3]);   
    writearray(out, sol.wdg, sol.nsize[4]);   
    
    // Close file:
    out.close();    
}

void readsolstruct(string filename, solstruct &sol, appstruct &app, masterstruct &master, string filemesh, Int mpirank)
{
    // Open file to read
    ifstream inmesh(filemesh.c_str(), ios::in | ios::binary);
    if (!inmesh) error("Unable to open file " + filemesh);    
    int *lsize = readiarrayfromdouble(inmesh, 1);
    int *nsize = readiarrayfromdouble(inmesh, lsize[0]);
    int *ndims = readiarrayfromdouble(inmesh, nsize[0]);
    Int ne = ndims[1];
    free(lsize);
    free(nsize);
    free(ndims);
    inmesh.close();            
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    //printf("Read sol struct from files...\n");  

    //dstype *tmp; // = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
    Int npe = master.ndims[5];
    Int npf = master.ndims[6];
    Int nc = app.ndims[5];
    Int ncu = app.ndims[6];
    Int nco = app.ndims[9];
    Int ncx = app.ndims[11];
    Int ncw = app.ndims[13];
    
    sol.lsize = readiarrayfromdouble(in, 1);
    sol.nsize = readiarrayfromdouble(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdouble(in, sol.nsize[0]);  
        
    //printf("%d %d %d %d %d %d %d %d %d\n", mpirank, sol.nsize[0], sol.nsize[1], sol.nsize[2], npe, ncx, nc, ncu, ne);
    
    if (sol.nsize[1] == npe*ncx*ne) {      
        if (mpirank==0) printf("Reading xdg from binary files \n");  
        readarray(in, &sol.xdg, sol.nsize[1]);
        sol.szxdg = sol.nsize[1];
    }
    else
        error("Input files are incorrect");
        
    if (sol.nsize[2] == npe*nc*ne) {
        if (mpirank==0) printf("Reading (u,q) from binary files \n");  
        readarray(in, &sol.udg, sol.nsize[2]);    
        sol.szudg = sol.nsize[2];
    }
    else if (sol.nsize[2] == npe*ncu*ne) {
        if (mpirank==0) printf("Reading u from binary files \n");  
        dstype *tmp;
        readarray(in, &tmp,     sol.nsize[2]);    
        sol.udg = (dstype*) malloc (sizeof (dstype)*npe*nc*ne);    
        cpuArraySetValue(sol.udg, zero, npe*nc*ne);

        cpuArrayInsert(sol.udg, tmp, npe, nc, ne, 0, npe, 0, ncu, 0, ne); 
        sol.nsize[2] = npe*nc*ne;     
        sol.szudg = sol.nsize[2];
        CPUFREE(tmp);
    }
    else if (sol.nsize[2] == 0) {
        if (mpirank==0) printf("compute the initial solution \n");  
        sol.udg = (dstype*) malloc (sizeof (dstype)*npe*nc*ne);    
        cpuArraySetValue(sol.udg, zero, npe*nc*ne);
        
        if (app.flag[1]==0) { //            
            cpuInituDriver(sol.udg, sol.xdg, app, ncx, nc, npe, ne, 0);    
        }
        else // wave problem
            cpuInitudgDriver(sol.udg, sol.xdg, app, ncx, nc, npe, ne, 0);                    
        sol.nsize[2] = npe*nc*ne;       
        sol.szudg = sol.nsize[2];
    }
    else
        error("Input files are incorrect");        
    
    #ifdef HAVE_ENZYME
        sol.dudg = (dstype*) malloc (sizeof (dstype)*npe*nc*ne);
        cpuArraySetValue(sol.dudg, zero, npe*nc*ne);
    #endif

    if ((sol.nsize[3] >0) && (sol.nsize[3] == npe*nco*ne)) {
        if (mpirank==0) printf("Reading vdg from binary files \n");  
        readarray(in, &sol.odg, sol.nsize[3]);    
        sol.szodg = sol.nsize[3];
    }
    else if (nco>0) {
        sol.odg = (dstype*) malloc (sizeof (dstype)*npe*nco*ne);
        cpuInitodgDriver(sol.odg, sol.xdg, app, ncx, nco, npe, ne, 0);       
        sol.nsize[3] = npe*nco*ne;
        sol.szodg = sol.nsize[3];
    } 
    #ifdef HAVE_ENZYME
        sol.dodg = (dstype*) malloc (sizeof (dstype)*npe*nco*ne);
        cpuArraySetValue(sol.dodg, zero, npe*nco*ne);
    #endif

    if ((sol.nsize[4] >0) && (sol.nsize[4] == npe*ncw*ne)) {
        if (mpirank==0) printf("Reading wdg from binary files \n");  
        readarray(in, &sol.wdg, sol.nsize[4]);    
        sol.szwdg = sol.nsize[4];
    }
    else if (ncw>0) {
        sol.wdg = (dstype*) malloc (sizeof (dstype)*npe*ncw*ne);    
        cpuInitwdgDriver(sol.wdg, sol.xdg, app, ncx, ncw, npe, ne, 0);        
        sol.nsize[4] = npe*ncw*ne;
        sol.szwdg = sol.nsize[4];
    }
    if (sol.nsize[5] >0) {
        if (mpirank==0) printf("Reading uh from binary files \n");  
        readarray(in, &sol.uh, sol.nsize[5]);
        app.read_uh = 1;
        sol.szuh = sol.nsize[5];
    }
    if (sol.nsize[6] > 0) {
        if (mpirank==0) printf("Reading xcg from binary files \n");  
        readarray(in, &sol.xcg, sol.nsize[6]);
        sol.szxcg = sol.nsize[6];
    }    
        
    in.close();            
}

void readInput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string filein, 
        Int mpiprocs, Int mpirank, Int fileoffset, Int omprank) 
{   
    if (mpirank==0) printf("Reading app from binary files \n");  
    string fileapp = filein + "app.bin";        
    readappstruct(fileapp,app);        
        
    if (mpirank==0) printf("Reading master from binary files \n");    
    string filemaster = filein + "master.bin";           
    readmasterstruct(filemaster, master);
    
    Int nd = master.ndims[0];    
    app.porder = (Int*) malloc(nd*sizeof(Int));
    app.szporder = nd;
    for (Int i=0; i<nd; i++)
        app.porder[i] = master.ndims[3];
    app.comm = (Int*) malloc(2*sizeof(Int));
    app.comm[0] = mpirank;
    app.comm[1] = mpiprocs;
    app.szcomm = 2;
                    
    // read meshsol structure
    if (mpiprocs>1) {     
        Int filenumber = mpirank+1-fileoffset; //file number     
        string filesol = filein + "sol" + NumberToString(filenumber) + ".bin";                    
        string filemesh = filein + "mesh" + NumberToString(filenumber) + ".bin";                    
        
        if (mpirank==0) printf("Reading initial solution from binary files \n");         
        readsolstruct(filesol, sol, app, master, filemesh, mpirank);    
        
        if (mpirank==0) printf("Reading mesh from binary files \n");            
        readmeshstruct(filemesh, mesh, sol, app, master, mpirank);                              
    }
    else {
        string filesol = filein + "sol.bin";              
        string filemesh = filein + "mesh.bin";                    

        if (mpirank==0) printf("Reading initial solution from binary files \n");                       
        readsolstruct(filesol, sol, app, master, filemesh, mpirank);    
        
        if (mpirank==0) printf("Reading mesh from binary files \n");         
        readmeshstruct(filemesh, mesh, sol, app, master, mpirank);                      
    }    
}

void writeOutput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string fileout, 
        Int mpiprocs, Int mpirank, Int fileoffset, Int omprank) 
{
    string fileoutapp = fileout + "app.bin";        
    writeappstruct(fileoutapp,app);
    
    string fileoutmaster = fileout + "master.bin";        
    writemasterstruct(fileoutmaster,master);
                
    // read meshsol structure
    if (mpiprocs>1) {     
        Int filenumber = mpirank+1-fileoffset; //file number
        
        string fileoutmesh = fileout + "mesh" + NumberToString(filenumber) + ".bin";                
        writemeshstruct(fileoutmesh,mesh);    
        
        string fileoutsol = fileout + "sol" + NumberToString(filenumber) + ".bin";                
        writesolstruct(fileoutsol,sol);    
    }
    else {
        string fileoutmesh = fileout + "mesh.bin";                
        writemeshstruct(fileoutmesh,mesh);    
        
        string fileoutsol = fileout + "sol.bin";                
        writesolstruct(fileoutsol,sol);    
    }
}

#endif    

