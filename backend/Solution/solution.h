#ifndef __SOLUTION_H__
#define __SOLUTION_H__

class CSolution {
private:
public:
    CDiscretization disc;  // spatial discretization class
    CPreconditioner prec;  // precondtioner class 
    CSolver solv;          // linear and nonlinear solvers
    //ofstream filestr;     // storing residual norms
    
    // constructor 
    CSolution(string filein, string fileout, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank, Int backend)   
       : disc(filein, fileout, mpiprocs, mpirank, fileoffset, omprank, backend),
         prec(disc, backend),
         solv(disc, backend) { 
       
//           disc.common.printinfo();
//           disc.app.printinfo();
//           disc.res.printinfo();
//           disc.tmp.printinfo();
//           disc.sol.printinfo();
//           disc.mesh.printinfo();
//           disc.master.printinfo();       
//           prec.precond.printinfo();
//           solv.sys.printinfo();
//           error("here");
       };        
    
    // destructor        
    ~CSolution(){  }; 

    // solve steady-state problems
    void SteadyProblem(Int backend);
    
    void SteadyProblem(ofstream &out, Int backend);    

    void SteadyProblem_PTC(ofstream &out, Int backend);    
            
    void SteadyProblem(CSolution *subprob, Int backend);    
    
    // advance the solution to the next time step using DIRK and BDF schemes
    void TimeStepping(ofstream &out, dstype time, Int istep, Int backend);
        
    // solve time-dependent problems
    void UnsteadyProblem(ofstream &out, Int backend);    
        
    // solve time-dependent problems
    void DIRK(Int backend);
    
    void DIRK(ofstream &out, Int backend);    
    
    void DIRK(CSolution *subprob, Int backend);

    // precompute some quantities
    void InitSolution(Int backend);    
    
    // solve both steady-state and time-dependent problems
    void SolveProblem(Int backend);    
    
    void SolveProblem(ofstream &out, Int backend);    
    
    void SolveProblem(CSolution *subprob, Int backend);    
    
    // save solutions in binary files
    void SaveSolutions(Int backend);    

    // save solutions in binary files
    void SaveSolutionsOnBoundary(Int backend);    
    
    // save solutions in binary files
    void SaveNodesOnBoundary(Int backend);    
    
    // read solutions from binary files
    void ReadSolutions(Int backend);        
    
    // save output in binary files
    void SaveOutputDG(Int backend);        
    
    // save output in binary files
    void SaveOutputCG(Int backend);        
};

#endif        


