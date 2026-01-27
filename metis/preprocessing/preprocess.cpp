// cmake -S . -B build
// cmake --build build -j4

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <array>
#include <filesystem>

using namespace std;

#define HAVE_MPI

#include <metis.h>
#include <parmetis.h>

#include "structs.hpp"
#include "TextParser.hpp"

#include "comparestructs.cpp"
#include "tinyexpr.cpp"
#include "helpers.cpp"
#include "readpdeapp.cpp"
#include "readmesh.cpp"
#include "makemesh.cpp"
#include "makemaster.cpp"
#include "domaindecomposition.cpp"
#include "connectivity.cpp"
#include "writebinaryfiles.cpp"
#include "parmetisexasim.cpp"

int main(int argc, char** argv)
{
    // ----------------------------------------------------
    // 0. Initialize MPI
    // ----------------------------------------------------
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc < 2) {
        if (rank==0) std::cerr << "Usage: ./parseinput <pdeapp.txt>\n";
        return 1;
    }    

    if (std::filesystem::exists(argv[1])) {
        if (rank==0) std::cout << "Generating Exasim's input files for this text file ("<< argv[1] << ") ... \n\n";
    } else {
        error("Error: Input file does not exist.\n");        
    }          
           
    InputParams params = parseInputFile(argv[1], rank);                           
    PDE pde = initializePDE(params, rank);    
     
    ParsedSpec spec = TextParser::parseFile(make_path(pde.datapath, pde.modelfile));        
    spec.exasimpath = pde.exasimpath;    

      // Mesh mesh = initializeMesh(params, pde);        
      // Master master = initializeMaster(pde, mesh);                                    
      // writeBinaryFiles(pde, mesh, master, spec);
  
    if (size == 1) {
      Mesh mesh = initializeMesh(params, pde);        
      Master master = initializeMaster(pde, mesh);                                    
      writeBinaryFiles(pde, mesh, master, spec);
    }
    else {
      Mesh mesh = initializeParMesh(params, spec, pde, MPI_COMM_WORLD);   

      Master master = initializeMaster(pde, mesh, rank);    
      
      //printf("%d %d %d %d\n", rank, pde.nsca, pde.nvec, pde.nten);

      if (rank==0) {
        writepde(pde, make_path(pde.datainpath, "app.bin"));
        writemaster(master, make_path(pde.datainpath, "master.bin"));    
      }
      MPI_Barrier(MPI_COMM_WORLD);

      callParMetis(mesh, pde, MPI_COMM_WORLD); 
      
      DMD dmd = initializeDMD(mesh, master, pde, MPI_COMM_WORLD);     

      writemesh(mesh, dmd, pde, master, MPI_COMM_WORLD);

      writesol(mesh, dmd, pde, master, MPI_COMM_WORLD);
      
      int nfe = mesh.nfe;
      int ne = dmd.elempart.size();    

      vector<double> xdg(master.npe*mesh.dim*ne, 0);
      select_columns(xdg.data(), mesh.xdg.data(), dmd.elempart_local.data(), master.npe*mesh.dim, mesh.ne);
      sendrecvdata(MPI_COMM_WORLD, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                   dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*mesh.dim);

      // if (pde.udgfile != "") {
      //   readParFieldFromBinaryFile(pde.udgfile, mesh.epart_local, mesh.udg, mesh.udgdims);
      //   int nc = mesh.udgdims[1];
      //   xdg.resize(master.npe*nc*ne, 0);
      //   select_columns(xdg.data(), mesh.udg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      //   sendrecvdata(MPI_COMM_WORLD, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
      //              dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      // }
      // if (pde.vdgfile != "") {
      //   readParFieldFromBinaryFile(pde.vdgfile, mesh.epart_local, mesh.vdg, mesh.vdgdims);
      //   int nc = mesh.vdgdims[1];
      //   xdg.resize(master.npe*nc*ne, 0);
      //   select_columns(xdg.data(), mesh.vdg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      //   sendrecvdata(MPI_COMM_WORLD, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
      //              dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      // }
      // if (pde.wdgfile != "") {
      //   readParFieldFromBinaryFile(pde.wdgfile, mesh.epart_local, mesh.wdg, mesh.wdgdims);
      //   int nc = mesh.wdgdims[1];
      //   xdg.resize(master.npe*nc*ne, 0);
      //   select_columns(xdg.data(), mesh.wdg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      //   sendrecvdata(MPI_COMM_WORLD, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
      //              dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      // }

      //------------------------------------------------------// 
      
      InputParams params1 = parseInputFile(argv[1], 0);                           
      PDE pde1 = initializePDE(params1, 0);         
      ParsedSpec spec1 = TextParser::parseFile(make_path(pde1.datapath, pde1.modelfile));        
      spec1.exasimpath = pde1.exasimpath;    
      Mesh mesh1 = initializeMesh(params1, pde1);        
      Master master1 = initializeMaster(pde1, mesh1);                                    
      vector<DMD> dmd1 = buildMeshDMD(pde1, mesh1, master1, spec1, rank); 

      comparePDE(pde1, pde, true, 1e-10);
      compareMaster(master1, master, true, 1e-10);        
      compareDMD(dmd1[rank], dmd, true);
      
      // std::vector<int> intelem;
      // intelem.reserve(mesh1.ne / size);
      // for (int e = 0; e < mesh1.ne; ++e) {
      //     if (mesh1.elem2cpu[e] == rank) intelem.push_back(e);
      // }
      // vector<int> t2ti(mesh.nfe*mesh.ne); 
      // select_columns(t2ti.data(), mesh1.t2t.data(), intelem.data(), mesh.nfe, mesh.ne);       
      // for (int i = 0; i < mesh.nfe*mesh.ne; i++) t2ti[i] = (t2ti[i] < 0) ? mesh.t2t[i] : t2ti[i];
      // //print2iarray(t2ti.data(), mesh.nfe, mesh.ne);
      // if (compareVecInt(mesh.t2t, t2ti, "t2t", true)) cout<<"Rank: "<<rank<<", t2t arrays are identical"<<endl;

      // vector<int> t2ti(mesh.nfe*mesh.ne); 
      // select_columns(t2ti.data(), mesh1.t2t.data(), dmd1[rank].elempart.data(), mesh.nfe, mesh.ne); 
      // vector<int> t2tj(mesh.nfe*mesh.ne); 
      // select_columns(t2tj.data(), mesh.t2t.data(), dmd.elempart_local.data(), mesh.nfe, mesh.ne);       
      // for (int i = 0; i < mesh.nfe*mesh.ne; i++) t2ti[i] = (t2ti[i] < 0) ? t2tj[i] : t2ti[i];      
      // if (compareVecInt(t2tj, t2ti, "t2t", true)) cout<<"Rank: "<<rank<<", t2t arrays are identical"<<endl;
      // 
      // t2ti.resize(nfe*ne);
      // t2tj.resize(nfe*ne);
      // 
      // // if (rank==0) {
      // //   printf("%d %d %d %d\n", rank, nfe, mesh.ne, ne);
      // //   print2iarray(dmd.nbsd.data(), 1, dmd.nbsd.size());
      // //   print2iarray(dmd.elemsendpts.data(), 1, dmd.elemsendpts.size());
      // //   print2iarray(dmd.elemrecvpts.data(), 1, dmd.elemrecvpts.size());
      // //   print2iarray(elemsend.data(), 1, nsend);
      // //   print2iarray(elemrecv.data(), 1, nrecv);
      // // }
      // sendrecvdata(MPI_COMM_WORLD, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
      //              dmd.localelemsend, dmd.localelemrecv, t2tj, t2tj, nfe);
      // 
      // select_columns(t2ti.data(), mesh1.t2t.data(), dmd1[rank].elempart.data(), nfe, ne); 
      // for (int i = 0; i < nfe*ne; i++) t2ti[i] = (t2ti[i] < 0) ? t2tj[i] : t2ti[i];      
      // if (compareVecInt(t2tj, t2ti, "t2t", true)) cout<<"Rank: "<<rank<<", resized t2t arrays are identical"<<endl;
      // 
      // int n0 = dmd1[rank].elempartpts[0];
      // int n1 = dmd1[rank].elempartpts[1];
      // int n2 = dmd1[rank].elempartpts[2];
      // vector<int> part0(n0, 0);
      // vector<int> part1(n1, 0);
      // vector<int> part2(n2, 0);
      // for (int i=0; i<n0; i++) part0[i] = dmd1[rank].elempart[i];
      // for (int i=0; i<n1; i++) part1[i] = dmd1[rank].elempart[n0+i];
      // for (int i=0; i<n2; i++) part2[i] = dmd1[rank].elempart[n0+n1+i];
      // 
      // int m0 = dmd.elempartpts[0];
      // int m1 = dmd.elempartpts[1];
      // int m2 = dmd.elempartpts[2];
      // vector<int> qart0(m0, 0);
      // vector<int> qart1(m1, 0);
      // vector<int> qart2(m2, 0);
      // for (int i=0; i<m0; i++) qart0[i] = dmd.elempart[i];
      // for (int i=0; i<m1; i++) qart1[i] = dmd.elempart[m0+i];
      // for (int i=0; i<m2; i++) qart2[i] = dmd.elempart[m0+m1+i];      
      // 
      // if (compareVecInt(qart0, part0, "interiorGlobal", true)) 
      //   cout<<"Rank: "<<rank<<", interiorGlobal arrays are identical"<<endl;
      // 
      // if (compareVecInt(qart1, part1, "interfaceGlobal", true)) 
      //   cout<<"Rank: "<<rank<<", interfaceGlobal arrays are identical"<<endl;
      // 
      // if (compareVecInt(qart2, part2, "exteriorGlobal", true)) 
      //   cout<<"Rank: "<<rank<<", exteriorGlobal arrays are identical"<<endl;
      
      if (compareVecInt(dmd.elempart, dmd1[rank].elempart, "elempart", true)) 
        cout<<"Rank: "<<rank<<", elempart arrays are identical"<<endl;

      if (compareVecInt(dmd.elempartpts, dmd1[rank].elempartpts, "elempartpts", true)) 
        cout<<"Rank: "<<rank<<", elempartpts arrays are identical"<<endl;

      if (compareVecInt(dmd.intepartpts, dmd1[rank].intepartpts, "intepartpts", true)) 
        cout<<"Rank: "<<rank<<", intepartpts arrays are identical"<<endl;

      // // if (rank==0) {
      // //   cout<<dmd.intepartpts.size()<<endl;
      // //   print2iarray(dmd.intepartpts.data(), 1, dmd.intepartpts.size());
      // //   print2iarray(dmd1[rank].intepartpts.data(), 1, dmd1[rank].intepartpts.size());
      // // }
      // 
      // for (int i=0; i<n0; i++) part0[i] = dmd1[rank].elem2cpu[i];
      // for (int i=0; i<n1; i++) part1[i] = dmd1[rank].elem2cpu[n0+i];
      // for (int i=0; i<n2; i++) part2[i] = dmd1[rank].elem2cpu[n0+n1+i];
      // for (int i=0; i<m0; i++) qart0[i] = dmd.elem2cpu[i];
      // for (int i=0; i<m1; i++) qart1[i] = dmd.elem2cpu[m0+i];
      // for (int i=0; i<m2; i++) qart2[i] = dmd.elem2cpu[m0+m1+i];      
      // 
      // if (compareVecInt(qart0, part0, "interior elem2cpu", true)) 
      //   cout<<"Rank: "<<rank<<", interior elem2cpu arrays are identical"<<endl;
      // 
      // if (compareVecInt(qart1, part1, "interface elem2cpu", true)) 
      //   cout<<"Rank: "<<rank<<", interface elem2cpu arrays are identical"<<endl;
      // 
      // if (compareVecInt(qart2, part2, "exterior elem2cpu", true)) 
      //   cout<<"Rank: "<<rank<<", exterior elem2cpu arrays are identical"<<endl;
      // 
      // if (compareArray3(dmd.elemrecv, dmd1[rank].elemrecv, "elemrecv", true))
      //   cout<<"Rank: "<<rank<<", elemrecv arrays are identical"<<endl;
      // 
      // if (compareArray3(dmd.elemsend, dmd1[rank].elemsend, "elemsend", true))
      //   cout<<"Rank: "<<rank<<", elemsend arrays are identical"<<endl;
      // else {
      //   if (rank==0) {
      //     printElemRecv(dmd.elemsend, "elemsend");
      //     printElemRecv(dmd1[rank].elemsend, "elemsend");        
      //   }
      // }
      // 
      // if (compareVecInt(dmd.elemrecvpts, dmd1[rank].elemrecvpts, "elemrecvpts", true)) 
      //   cout<<"Rank: "<<rank<<", elemrecvpts arrays are identical"<<endl;
      // 
      // if (compareVecInt(dmd.elemsendpts, dmd1[rank].elemsendpts, "elemsendpts", true)) 
      //   cout<<"Rank: "<<rank<<", elemsendpts arrays are identical"<<endl;
      // else {
      //   if (rank==0) {
      //     print2iarray(dmd.elemsendpts.data(), 1, dmd.elemsendpts.size());
      //     print2iarray(dmd1[rank].elemsendpts.data(), 1, dmd1[rank].elemsendpts.size());
      //   }
      // }
      // 
      // if (compareVecInt(dmd.nbsd, dmd1[rank].nbsd, "nbsd", true)) 
      //   cout<<"Rank: "<<rank<<", nbsd arrays are identical"<<endl;
      // 
      // // vector<int> tm(mesh.nfe*mesh.ne, 0); 
      // // for (int i = 0; i < mesh.nfe*mesh.ne; i++) tm[i] = dmd1[rank].bf[i];
      // // if (compareVecInt(mesh.bf, tm, "bf", true)) 
      // //   cout<<"Rank: "<<rank<<", bf arrays are identical"<<endl;
      if (compareVecInt(mesh.bf, dmd1[rank].bf, "bf", true)) 
        cout<<"Rank: "<<rank<<", bf arrays are identical"<<endl;

      // if (rank==0) {
      //   print2iarray(mesh.bf.data(), mesh.nfe, dmd.elempart.size());   
      //   print2iarray(dmd1[rank].bf.data(), mesh.nfe, dmd.elempart.size());   
      // }
      
      // 
      // vector<int> ti(mesh.nve*ne); 
      // select_columns(ti.data(), mesh1.t.data(), dmd1[rank].elempart.data(), mesh.nve, ne); 
      // if (compareVecInt(mesh.tg, ti, "tglobal", true)) 
      //   cout<<"Rank: "<<rank<<", tglobal arrays are identical"<<endl;
      
      vector<double> xdg1(master.npe*mesh.dim*ne, 0);
      select_columns(xdg1.data(), mesh1.xdg.data(), dmd1[rank].elempart.data(), master.npe*mesh.dim, ne);
      if (compareVecDouble(xdg, xdg1, "xdg", true, 1e-7)) 
        cout<<"Rank: "<<rank<<", xdg arrays are identical"<<endl;
       
      // if (rank==0) {        
      //   //print2darray(xdg.data(), master.npe*mesh.dim, mesh.ne);        
      //   //print2darray(xdg1.data(), master.npe*mesh.dim, mesh.ne);        
      //   //print2iarray(dmd.elempart_local.data(), 1, mesh.ne);        
      //   //print2iarray(mesh.tg.data(), mesh.nve, mesh.ne);        
      //   //print2iarray(ti.data(), mesh.nve, mesh.ne);        
      //   //print2iarray(mesh.t2t.data(), mesh.nfe, mesh.ne);        
      //   //print2iarray(t2ti.data(), mesh.nfe, mesh.ne);        
      //   //print2iarray(mesh.bf.data(), mesh.nfe, mesh.ne);        
      //   //print2iarray(tm.data(), mesh.nfe, mesh.ne);
      //   // vector<int> fi(mesh.nfe*mesh.ne);      
      //   // select_columns(fi.data(), mesh1.f.data(), dmd1[rank].elempart.data(), mesh.nfe, mesh.ne);       
      //   // print2iarray(fi.data(), mesh.nfe, mesh.ne);
      // 
      //   //printElemRecv(dmd.elemrecv);
      //   //printElemRecv(dmd1[rank].elemrecv);
      //   //printElemRecv(dmd.elemsend, "elemsend");
      //   //printElemRecv(dmd1[rank].elemsend, "elemsend");
      //   //print2iarray(dmd.elempartpts.data(), 1, dmd.elempartpts.size());        
      // }
      // 
      // // if (rank==0) {
      // //   print2iarray(part1.data(), 1, part1.size());  
      // //   print2iarray(elemclass.interfaceGlobal.data(), 1, elemclass.interfaceGlobal.size());  
      // //   print2iarray(part2.data(), 1, part2.size());  
      // //   print2iarray(elemclass.neighborElemGlobal.data(), 1, elemclass.neighborElemGlobal.size());  
      // // }
      // 
      // //print2iarray(dmd1[rank].elempartpts.data(), 1, dmd1[rank].elempartpts.size());  
      // 
      // // print2iarray(mesh.t2t.data(), mesh.nfe, mesh.ne);
      // // vector<int> t2ti(mesh.nfe*mesh.ne); 
      // // select_columns(t2ti.data(), mesh1.t2t.data(), dmd1[rank].elempart.data(), mesh.nfe, mesh.ne); 
      // // print2iarray(t2ti.data(), mesh.nfe, mesh.ne);
      // 
      // // compareDMD(dmd1[rank], dmd, true);
      // // print2iarray(dmd1[rank].elempart.data(), 1, dmd1[rank].elempart.size());
      // // print2iarray(dmd.elempart.data(), 1, dmd.elempart.size());
      // //print2iarray(mesh1.elem2cpu.data(), 1, mesh1.elem2cpu.size());
      // 
      //  //printvector(mesh.epart_local, "mesh.epart_local", MPI_COMM_WORLD);

    }
    
    if (rank==0) std::cout << "\n******** Done with generating input files for EXASIM ********\n";
  
    MPI_Finalize();
    return 0;
}

