#ifndef __ERRORMSG
#define __ERRORMSG

static void PrintErrorAndExit(const char* errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg, file, line );
    
#ifdef  HAVE_MPI       
    MPI_Finalize();    
#endif
    
    exit( 1 );    
}

static void PrintErrorAndExit(string errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg.c_str(), file, line );
    
#ifdef  HAVE_MPI       
    MPI_Finalize();    
#endif
    
    exit( 1 );    
}

static void PrintMsg(const char* errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg, file, line );    
}

static void PrintMsg(string errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg.c_str(), file, line );
}

#define error( errmsg ) (PrintErrorAndExit( errmsg, __FILE__, __LINE__ ))
#define display( errmsg ) (PrintMsg( errmsg, __FILE__, __LINE__ ))

#endif
