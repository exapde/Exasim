#ifndef __READBINARYFILES
#define __READBINARYFILES

struct appsolstruct {     
    string exasimpath = "";  
    string filein = "";       // Name of binary file with input data
    string fileout = "";      // Name of binary file to write the solution                
    int modelnumber;
    int mpiRank = 0;
    int mpiProcs = 0;
    int saveSolFreq = 0;
    int saveSolOpt = 0;
    int saveSolBouFreq = 0;
    int tdep = 0;
    int wave = 0;
    int currentstep = 0;
    int tsteps = 0;
    int timestepOffset = 0;
  
    int nd = 0;
    int nc = 0;
    int ncu = 0;
    int ncq = 0;
    int ncw = 0;
    int nce = 0;
    int ncx = 0;
    int nco = 0;
    int szxcg = 0;    
    int nsca = 0;
    int nvec = 0;
    int nten = 0;
    int nsurf = 0;
    int nvqoi = 0;    
    int ne1 = 0;
    int elemtype = 0;     
    int nodetype = 0;     
    int porder = 0;
    int pgauss = 0;
    int npe = 0;
    int npf = 0;
    int nge = 0;
    int ngf = 0;

    int* flag=nullptr;
    int* problem=nullptr;
    
    int *facecon=nullptr;    // face-to-element connectivities     
    int *elemcon=nullptr;    // element-to-face connectivities
    int *perm=nullptr;       // indices of element nodes on faces
    int *eblks=nullptr;    // element blocks
    int *fblks=nullptr;    // face blocks    
    int *nbsd=nullptr;
    int *elemsend=nullptr;
    int *elemrecv=nullptr;
    int *elemsendpts=nullptr;
    int *elemrecvpts=nullptr;
    int *elempart=nullptr;
    int *elempartpts=nullptr;
    int *cgelcon=nullptr;
    int *rowent2elem=nullptr;
    int *cgent2dgent=nullptr;
    int *colent2elem=nullptr;  
    int *rowe2f1=nullptr;
    int *cole2f1=nullptr;
    int *ent2ind1=nullptr;
    int *rowe2f2=nullptr;
    int *cole2f2=nullptr;
    int *ent2ind2=nullptr;
    int *perm=nullptr;
        
    dstype time;
    dstype* dt=nullptr;    
    dstype* factor=nullptr;
    dstype* physicsparam=nullptr;
    dstype* externalparam=nullptr;

    dstype *xdg=nullptr; // spatial coordinates
    dstype *udg=nullptr; // solution (u, q) 
    dstype *odg=nullptr; // auxilary term 
    dstype *wdg=nullptr; // dw/dt = u (wave problem)  
};

void readappstruct(string filename, appsolstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) error("Unable to open file " + filename);
           
    int* lsize=nullptr;
    int* nsize=nullptr;
    int* ndims=nullptr;
  
    /* Read data to app structure */            
    lsize = readiarrayfromdouble(in, 1);
    nsize = readiarrayfromdouble(in, lsize[0]);
    ndims = readiarrayfromdouble(in, nsize[0]);
    app.flag = readiarrayfromdouble(in, nsize[1]);
    app.problem = readiarrayfromdouble(in, size[2]);
    readarray(in, &app.externalparam, nsize[3]);
    readarray(in, &app.dt, nsize[4]);                
    readarray(in, &app.factor, nsize[5]);       
    readarray(in, &app.physicsparam, nsize[6]);       
    
    app.nc = ndims[5]; // number of compoments of (u, q)
    app.ncu = ndims[6];// number of compoments of (u)        
    app.ncq = ndims[7];// number of compoments of (q)
    app.nco = ndims[9];// number of compoments of (o)    
    app.ncx = ndims[11];// number of compoments of (xdg)        
    app.nce = ndims[12];// number of compoments of (output)        
    app.ncw = ndims[13];//number of compoments of (w)
    app.nsca = ndims[14];// number of components of scalar fields for visualization
    app.nvec = ndims[15];// number of components of vector fields for visualization
    app.nten = ndims[16];// number of components of tensor fields for visualization
    app.nsurf = ndims[17];// number of components of surface fields for visualization, storage, and QoIs
    app.nvqoi = ndims[18];// number of volume quantities of interest (QoIs)    
  
    app.tdep = app.flag[0];      // 0: steady-state; 1: time-dependent;  
    app.wave = app.flag[1];
  
    app.saveSolFreq = app.problem[17];    
    app.saveSolOpt = app.problem[18];    
    app.timestepOffset = app.problem[19];    
    app.saveSolBouFreq = app.problem[21];   
  
    app.tsteps = nsize[4];  // number of time steps          
    
    free(lsize);
    free(nsize);
    free(ndims);

    // Close file:
    in.close();
}

void readmasterstruct(string filename, appsolstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) error("Unable to open file " + filename);
        
    int* lsize=nullptr;
    int* nsize=nullptr;
    int* ndims=nullptr;
  
    /* Read data to app structure */        
    lsize = readiarrayfromdouble(in, 1);
    nsize = readiarrayfromdouble(in, lsize[0]);        
    ndims = readiarrayfromdouble(in, nsize[0]);

    app.nd = ndims[0];     // spatial dimension    
    app.elemtype = ndims[1]; 
    app.nodetype = ndims[2]; 
    app.porder = ndims[3]; 
    app.pgauss = ndims[4]; 
    app.npe = ndims[5]; // number of nodes on master element
    app.npf = ndims[6]; // number of nodes on master face       
    app.nge = ndims[7]; // number of gauss points on master element
    app.ngf = ndims[8]; // number of gauss points on master face          
  
    free(lsize);
    free(nsize);
    free(ndims);
  
    // Close file:
    in.close();
}

void readmeshstruct(string filename, appsolstruct &app)
{
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);

    if (!in) 
        error("Unable to open file " + filename);
        
    int* lsize=nullptr;
    int* nsize=nullptr;
    int* ndims=nullptr;

    lsize = readiarrayfromdouble(in, 1);
    nsize = readiarrayfromdouble(in, lsize[0]);
    ndims = readiarrayfromdouble(in, nsize[0]);
    app.facecon = readiarrayfromdouble(in, nsize[1]); //  
    app.eblks = readiarrayfromdouble(in, nsize[2]);  //
    app.fblks = readiarrayfromdouble(in, nsize[3]);  //
    app.nbsd = readiarrayfromdouble(in, nsize[4]);
    app.elemsend = readiarrayfromdouble(in, nsize[5]);
    app.elemrecv = readiarrayfromdouble(in, nsize[6]);
    app.elemsendpts = readiarrayfromdouble(in, nsize[7]);
    app.elemrecvpts = readiarrayfromdouble(in, nsize[8]);
    app.elempart = readiarrayfromdouble(in, nsize[9]);
    app.elempartpts = readiarrayfromdouble(in, nsize[10]);    
    app.cgelcon = readiarrayfromdouble(in, nsize[11]); //
    app.rowent2elem = readiarrayfromdouble(in, nsize[12]); //
    app.cgent2dgent = readiarrayfromdouble(in, nsize[13]); //
    app.colent2elem = readiarrayfromdouble(in, nsize[14]); //
    app.rowe2f1 = readiarrayfromdouble(in, app.nsize[15]);     //
    app.cole2f1 = readiarrayfromdouble(in, app.nsize[16]);    //
    app.ent2ind1 = readiarrayfromdouble(in, app.nsize[17]);   //
    app.rowe2f2 = readiarrayfromdouble(in, app.nsize[18]);    //
    app.cole2f2 = readiarrayfromdouble(in, app.nsize[19]);    //
    app.ent2ind2 = readiarrayfromdouble(in, app.nsize[20]);   //
    app.f2e = readiarrayfromdouble(in, app.nsize[21]);        //
    app.elemcon = readiarrayfromdouble(in, app.nsize[22]);    //
    app.perm = readiarrayfromdouble(in, app.nsize[23]);
    app.bf = readiarrayfromdouble(in, app.nsize[24]);
  
    // Close file:
    in.close();            

    free(lsize);
    free(nsize);
    free(ndims);  
}

void readInput(appstruct &app, masterstruct &master, meshstruct &mesh, string filein, 
        int mpiprocs, int mpirank, int fileoffset) 
{   
    if (mpirank==0) printf("Reading app from binary files \n");  
    string fileapp = filein + "app.bin";        
    readappstruct(fileapp,app);        
        
    if (mpirank==0) printf("Reading master from binary files \n");    
    string filemaster = filein + "master.bin";           
    readmasterstruct(filemaster, master);
    
    int nd = ndims[0];    
    app.porder = (int*) malloc(nd*sizeof(int));
    app.szporder = nd;
    for (int i=0; i<nd; i++)
        app.porder[i] = ndims[3];
    app.comm = (int*) malloc(2*sizeof(int));
    app.comm[0] = mpirank;
    app.comm[1] = mpiprocs;
    app.szcomm = 2;
                    
    // read meshsol structure
    if (mpiprocs>1) {     
        int filenumber = mpirank+1-fileoffset; //file number     
        string filemesh = filein + "mesh" + NumberToString(filenumber) + ".bin";                    
                
        if (mpirank==0) printf("Reading mesh from binary files \n");            
        readmeshstruct(filemesh, mesh);                              
    }
    else {
        string filemesh = filein + "app.bin";                    
        
        if (mpirank==0) printf("Reading mesh from binary files \n");         
        readmeshstruct(filemesh, mesh); 
    }    
}


#endif    

