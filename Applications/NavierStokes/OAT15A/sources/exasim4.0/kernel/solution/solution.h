#ifndef __SOLUTION_H__
#define __SOLUTION_H__

class CSolution {
private:
public:
    CDiscretization disc;  // spatial discretization class
    CPreconditioner prec;  // precondtioner class 
    CSolver solv;          // linear and nonlinear solvers

    // constructor 
    CSolution(string filein, string fileout, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank, Int backend)   
       : disc(filein, fileout, mpiprocs, mpirank, ompthreads, omprank, backend),
         prec(disc, backend),
         solv(disc, backend){ };        
    
    // destructor        
    ~CSolution(){ }; 

    // solve steady-state problems
    void SteadyProblem(Int backend);    
            
    // advance the solution to the next time step using DIRK and BDF schemes
    void TimeStepping(dstype time, Int istep, Int backend);
        
    // solve time-dependent problems
    void UnsteadyProblem(Int backend);    
        
    // solve both steady-state and time-dependent problems
    void SolveProblem(Int backend);    
};

#endif        

