#ifndef __READBINARYFILES
#define __READBINARYFILES

#include "../AppDriver/inituDriver.cpp"
#include "../AppDriver/initqDriver.cpp"
#include "../AppDriver/initudgDriver.cpp"
#include "../AppDriver/initwdgDriver.cpp"
#include "../AppDriver/initodgDriver.cpp"

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

//     app.nproc = app.ndims[0];    // number of MPI ranks
//     app.nd = app.ndims[1]; // spatial dimension    
//     app.ne = app.ndims[2]; // total number of elements
//     app.nf = app.ndims[3]; // total number of faces
//     app.nv = app.ndims[4]; // total number of vertices      
//     app.nc = app.ndims[5]; // number of compoments of (u, q, p)
//     app.ncu = app.ndims[6];// number of compoments of (u)
//     app.ncq = app.ndims[7];// number of compoments of (q)
//     app.ncp = app.ndims[8];// number of compoments of (p)
//     app.nco = app.ndims[9];// number of compoments of (o)
//     app.nch = app.ndims[10];// number of compoments of (uhat)
//     app.ncx = app.ndims[11];// number of compoments of (xdg)
    
    Int i, ncu, ncq, ncp;
    ncu = app.ndims[6];// number of compoments of (u)
    ncq = app.ndims[7];// number of compoments of (q)
    ncp = app.ndims[8];// number of compoments of (p)
    
    if (ncu>0) {
        app.fc_u = (dstype*) malloc(sizeof(dstype)*ncu);
        app.dtcoef_u = (dstype*) malloc(sizeof(dstype)*ncu);
        for (i=0; i<ncu; i++) {
            app.fc_u[i] = 1.0;     //app.factor[i];
            app.dtcoef_u[i] = 1.0; //app.factor[i];
        }
    }        
    if (ncq>0) {
        app.fc_q = (dstype*) malloc(sizeof(dstype)*ncq);
        app.dtcoef_q = (dstype*) malloc(sizeof(dstype)*ncq);
        for (i=0; i<ncq; i++) {
            app.fc_q[i] = 1.0;     //app.factor[ncu+i];
            app.dtcoef_q[i] = 1.0; //app.factor[ncu+i];
        }
    }        
    if (ncp>0) {
        app.fc_p = (dstype*) malloc(sizeof(dstype)*ncp);
        app.fc_p = (dstype*) malloc(sizeof(dstype)*ncp);
        for (i=0; i<ncp; i++) {
            app.fc_p[i] = 1.0;    // app.factor[ncu+ncq+i];
            app.dtcoef_p[i] =1.0; // app.factor[ncu+ncq+i];
        }
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


void readmeshstruct(string filename, meshstruct &mesh)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    //printf("Read mesh struct from files...\n");  
    
    mesh.lsize = readiarrayfromdouble(in, 1);
    mesh.nsize = readiarrayfromdouble(in, mesh.lsize[0]);
    mesh.ndims = readiarrayfromdouble(in, mesh.nsize[0]);
    mesh.facecon = readiarrayfromdouble(in, mesh.nsize[1]);
    mesh.eblks = readiarrayfromdouble(in, mesh.nsize[2]);
    mesh.fblks = readiarrayfromdouble(in, mesh.nsize[3]);
    mesh.nbsd = readiarrayfromdouble(in, mesh.nsize[4]);
    mesh.elemsend = readiarrayfromdouble(in, mesh.nsize[5]);
    mesh.elemrecv = readiarrayfromdouble(in, mesh.nsize[6]);
    mesh.elemsendpts = readiarrayfromdouble(in, mesh.nsize[7]);
    mesh.elemrecvpts = readiarrayfromdouble(in, mesh.nsize[8]);
    mesh.elempart = readiarrayfromdouble(in, mesh.nsize[9]);
    mesh.elempartpts = readiarrayfromdouble(in, mesh.nsize[10]);    
    mesh.cgelcon = readiarrayfromdouble(in, mesh.nsize[11]);
    mesh.rowent2elem = readiarrayfromdouble(in, mesh.nsize[12]);
    mesh.cgent2dgent = readiarrayfromdouble(in, mesh.nsize[13]);
    mesh.colent2elem = readiarrayfromdouble(in, mesh.nsize[14]);
    mesh.rowe2f1 = readiarrayfromdouble(in, mesh.nsize[15]);
    mesh.cole2f1 = readiarrayfromdouble(in, mesh.nsize[16]);
    mesh.ent2ind1 = readiarrayfromdouble(in, mesh.nsize[17]);
    mesh.rowe2f2 = readiarrayfromdouble(in, mesh.nsize[18]);
    mesh.cole2f2 = readiarrayfromdouble(in, mesh.nsize[19]);
    mesh.ent2ind2 = readiarrayfromdouble(in, mesh.nsize[20]);
    
//     mesh.nd = mesh.ndims[0]; // spatial dimension    
//     mesh.ne = mesh.ndims[1]; // number of elements in this subdomain 
//     mesh.nf = mesh.ndims[2]; // number of faces in this subdomain 
//     mesh.nv = mesh.ndims[3]; // number of vertices in this subdomain       
//     mesh.nfe = mesh.ndims[4]; // number of faces per element        
//     mesh.nbe = mesh.ndims[5]; // number of blocks for elements 
//     mesh.neb = mesh.ndims[6]; // maximum number of elements per block
//     mesh.nbf = mesh.ndims[7]; // number of blocks for faces   
//     mesh.nfb = mesh.ndims[8]; // maximum number of faces per block                                        
    
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
    
//     sol.ne = sol.ndims[0]; // number of elements in this subdomain 
//     sol.nf = sol.ndims[1]; // number of faces in this subdomain 
//     sol.nfe = sol.ndims[2]; // number of faces per element        
//     sol.npe = sol.ndims[3];  // number of nodes on master element
//     sol.npf = sol.ndims[4];  // number of nodes on master face       
//     sol.nc = sol.ndims[5]; // number of compoments of (u, q, p)
//     sol.ncu = sol.ndims[6];// number of compoments of (u)
//     sol.ncq = sol.ndims[7];// number of compoments of (q)
//     sol.ncp = sol.ndims[8];// number of compoments of (p)
//     sol.nco = sol.ndims[9];// number of compoments of (o)
//     sol.nch = sol.ndims[10];// number of compoments of (uhat)
//     sol.ncx = sol.ndims[11];// number of compoments of (xdg)
    
    // Close file:
    in.close();            
}

void readsolstructManaged(string filename, solstruct &sol)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    printf("Read sol struct from files...\n");  

#ifndef HAVE_CUDA    
    sol.lsize = readiarrayfromdouble(in, 1);
    sol.nsize = readiarrayfromdouble(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdouble(in, sol.nsize[0]);    
    readarray(in, &sol.xdg, sol.nsize[1]);
    readarray(in, &sol.udg, sol.nsize[2]);    
    readarray(in, &sol.odg, sol.nsize[3]);    
    readarray(in, &sol.wdg, sol.nsize[4]);    
#endif
    
#ifdef HAVE_CUDA
    sol.lsize = readiarrayfromdoubleManaged(in, 1);
    sol.nsize = readiarrayfromdoubleManaged(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdoubleManaged(in, sol.nsize[0]);    
    readarrayManaged(in, &sol.xdg, sol.nsize[1]);
    readarrayManaged(in, &sol.udg, sol.nsize[2]);    
    readarrayManaged(in, &sol.odg, sol.nsize[3]);    
    readarrayManaged(in, &sol.wdg, sol.nsize[4]);    
#endif    
        
    // Close file:
    in.close();            
}

void readsolstructZeroCopy(string filename, solstruct &sol)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    printf("Read sol struct from files...\n");  

#ifndef HAVE_CUDA    
    sol.lsize = readiarrayfromdouble(in, 1);
    sol.nsize = readiarrayfromdouble(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdouble(in, sol.nsize[0]);    
    readarray(in, &sol.xdg, sol.nsize[1]);
    readarray(in, &sol.udg, sol.nsize[2]);    
    readarray(in, &sol.odg, sol.nsize[3]);    
    readarray(in, &sol.wdg, sol.nsize[4]);    
#endif
    
#ifdef HAVE_CUDA
    sol.lsize = readiarrayfromdoubleZeroCopy(in, 1);
    sol.nsize = readiarrayfromdoubleZeroCopy(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdoubleZeroCopy(in, sol.nsize[0]);    
    readarrayZeroCopy(in, &sol.xdg, sol.nsize[1]);
    readarrayZeroCopy(in, &sol.udg, sol.nsize[2]);    
    readarrayZeroCopy(in, &sol.odg, sol.nsize[3]);   
    readarrayZeroCopy(in, &sol.wdg, sol.nsize[4]);   
#endif    
        
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


void readsolstruct(string filename, solstruct &sol, appstruct &app, masterstruct &master, meshstruct &mesh)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
    
    //printf("Read sol struct from files...\n");  

    //dstype *tmp; // = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
    Int npe = master.ndims[5];
    Int nc = app.ndims[5];
    Int ncu = app.ndims[6];
    Int nco = app.ndims[9];
    Int ncx = app.ndims[11];
    Int ncw = app.ndims[13];
    Int ne = mesh.ndims[1];    

    sol.lsize = readiarrayfromdouble(in, 1);
    sol.nsize = readiarrayfromdouble(in, sol.lsize[0]);
    sol.ndims = readiarrayfromdouble(in, sol.nsize[0]);  
        
    sol.nsize[1] = npe*ncx*ne;
    sol.nsize[2] = npe*nc*ne;
    sol.nsize[3] = npe*nco*ne;
    sol.nsize[4] = npe*ncw*ne;
    
    readarray(in, &sol.xdg, sol.nsize[1]);
//    readarray(in, &sol.udg, sol.nsize[2]);
//     readarray(in, &tmp,     sol.nsize[2]);    
//     readarray(in, &sol.odg, sol.nsize[3]);  
//     readarray(in, &sol.wdg, sol.nsize[4]);              
            
    if (nc>0) {
        sol.udg = (dstype*) malloc (sizeof (dstype)*npe*nc*ne);    
        ArraySetValue(sol.udg, zero, npe*nc*ne, 0);
        if (app.flag[1]==0) {            
            InituDriver(sol.udg, sol.xdg, app, ncx, nc, npe, ne, 0);    
        }
        else 
            InitudgDriver(sol.udg, sol.xdg, app, ncx, nc, npe, ne, 0);            
    }
    if (nco>0) {
        sol.odg = (dstype*) malloc (sizeof (dstype)*npe*nco*ne);    
        InitwdgDriver(sol.odg, sol.xdg, app, ncx, nco, npe, ne, 0);        
    }
    if (ncw>0) {
        sol.wdg = (dstype*) malloc (sizeof (dstype)*npe*ncw*ne);    
        InitodgDriver(sol.wdg, sol.xdg, app, ncx, nco, npe, ne, 0);        
    }
            
//     // insert u into udg
//     ArrayInsert(sol.udg, tmp, npe, nc, ne, 0, npe, 0, ncu, 0, ne, 0); 
//     sol.nsize[2] = npe*nc*ne;
    
    //CPUFREE(tmp);
    
    // Close file:
    in.close();            
}

void readInput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string filein, 
        Int mpiprocs, Int mpirank, Int ompthreads, Int omprank) 
{   
    string fileapp = filein + "app.bin";        
    readappstruct(fileapp,app);        
        
    string filemaster = filein + "master.bin";           
    readmasterstruct(filemaster, master);
    
    Int nd = master.ndims[0];    
    app.porder = (Int*) malloc(nd*sizeof(Int));
    for (Int i=0; i<nd; i++)
        app.porder[i] = master.ndims[3];
    app.comm = (Int*) malloc(2*sizeof(Int));
    app.comm[0] = mpirank;
    app.comm[1] = mpiprocs;
                    
    // read meshsol structure
    if (mpiprocs>1) {     
        Int filenumber = mpirank+1; //file number     
        
        string filemesh = filein + "mesh" + NumberToString(filenumber) + ".bin";                    
        readmeshstruct(filemesh, mesh);              
                
        string filesol = filein + "sol" + NumberToString(filenumber) + ".bin";                    
        //readsolstruct(filesol, sol);    
        readsolstruct(filesol, sol, app, master, mesh);    
    }
    else {
        string filemesh = filein + "mesh.bin";                    
        readmeshstruct(filemesh, mesh);              
        
        string filesol = filein + "sol.bin";                    
        //readsolstruct(filesol, sol);              
        readsolstruct(filesol, sol, app, master, mesh);    
    }    
}

void writeOutput(appstruct &app, masterstruct &master, meshstruct &mesh, solstruct &sol, string fileout, 
        Int mpiprocs, Int mpirank, Int ompthreads, Int omprank) 
{
    string fileoutapp = fileout + "app.bin";        
    writeappstruct(fileoutapp,app);
    
    string fileoutmaster = fileout + "master.bin";        
    writemasterstruct(fileoutmaster,master);
                
    // read meshsol structure
    if (mpiprocs>1) {     
        Int filenumber = mpirank+1; //file number
        
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

