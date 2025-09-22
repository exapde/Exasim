#ifndef __ERRORMSG
#define __ERRORMSG

void error(const char* errstr)
{
#ifdef  HAVE_MPI            // TODO: Print only with master processor
    cout << "Error: ";
    cout << errstr;
    cout << endl;

//    free(sys.requests); free(sys.statuses);       // TODO: Free variables allocated dynamically
//    free(sys.requestsEnt2entWeight); free(sys.statusesEnt2entWeight);
//     MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#else
    cout << "Error: ";
    cout << errstr;
    cout << endl;
#endif

    exit(1);
}

void error(string errstr)
{
  error(errstr.c_str());
}

#endif
