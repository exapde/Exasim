/**
 * @file errormsg.cpp
 * @brief Provides error handling utilities for printing error messages and terminating the program.
 *
 * Defines functions to print error messages along with file name and line number,
 * and then exit the program. If MPI is enabled (HAVE_MPI defined), MPI_Finalize is called before exiting.
 *
 * Functions:
 * - PrintErrorAndExit(const char* errmsg, const char* file, int line): Prints a C-style error message and exits.
 * - PrintErrorAndExit(std::string errmsg, const char* file, int line): Prints a std::string error message and exits.
 *
 * Macros:
 * - error(errmsg): Calls PrintErrorAndExit with the provided error message, current file, and line.
 *
 * Usage:
 * Call error("Your error message") to print the message and terminate the program.
 *
 * @note The PrintMsg functions are commented out and not currently in use.
 */
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

// static void PrintMsg(const char* errmsg, const char *file, int line ) 
// {    
//     printf( "%s in %s at line %d\n", errmsg, file, line );    
// }
// 
// static void PrintMsg(string errmsg, const char *file, int line ) 
// {    
//     printf( "%s in %s at line %d\n", errmsg.c_str(), file, line );
// }

#define error( errmsg ) (PrintErrorAndExit( errmsg, __FILE__, __LINE__ ))
//#define display( errmsg ) (PrintMsg( errmsg, __FILE__, __LINE__ ))

#endif
