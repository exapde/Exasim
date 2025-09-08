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


std::string trim_dir(const std::string& s) {
    return std::filesystem::path{s}.parent_path().string();   // use .native() if you want OS-preferred slashes
}

bool ensure_dir(const std::string& dir) {
    std::filesystem::path p(dir);
    if (std::filesystem::exists(p)) return std::filesystem::is_directory(p);  // false if it's a file
    return std::filesystem::create_directories(p);               // creates parents as needed
}


std::string make_path(const std::string& str1, const std::string& str2) {
    std::filesystem::path base = str1;
    std::filesystem::path tail = str2;

    // If tail is absolute, strip its root so it becomes relative
    if (tail.is_absolute())
        tail = tail.relative_path();

    std::filesystem::path full = base / tail;
    return full.lexically_normal().string();
}

std::string trimToSubstringAtFirstOccurence(const std::string& fullPath, const std::string& keyword) {
    std::size_t pos = fullPath.find(keyword);  // Use find to get the first occurrence
    if (pos != std::string::npos) {
        return fullPath.substr(0, pos + keyword.length());
    }
    else {      
      return "";
    }
}

std::string trimToSubstringAtLastOccurence(const std::string& fullPath, const std::string& keyword) {
    std::size_t pos = fullPath.rfind(keyword);  // Use rfind to get the last occurrence
    if (pos != std::string::npos) {
        return fullPath.substr(0, pos + keyword.length());
    }
    else {      
      return "";
    }
}

#endif
