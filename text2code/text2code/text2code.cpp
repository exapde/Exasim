#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <array>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <cmath>

#ifdef HAVE_METIS
#include <metis.h>
#endif

// g++ -O2 -std=c++17 text2bina.cpp -o text2bina -lblas -llapack
// g++ -O2 -DHAVE_METIS -std=c++17 text2bina.cpp -o text2bina -lblas -llapack -I../METIS/build/xinclude -L../METIS/build/libmetis -lmetis -I../GKlib/include -L../GKlib/lib -lGKlib
// clang++ -fsanitize=address -fno-omit-frame-pointer -g -O2 -DHAVE_METIS -std=c++17 text2bina.cpp -o text2bina -lblas -llapack -I../METIS/build/xinclude -L../METIS/build/libmetis -lmetis -I../GKlib/include -L../GKlib/lib -lGKlib
// g++ -O2 -DHAVE_METIS -std=c++17 text2code.cpp -o text2code -lblas -llapack -I../METIS/build/xinclude -L../METIS/build/libmetis -lmetis -I../GKlib/include -L../GKlib/lib -lGKlib
// g++ -O2 -std=c++17 text2code.cpp -o text2code -lblas -llapack 

using namespace std;
        
#include "tinyexpr.cpp"
#include "helpers.cpp"
#include "readpdeapp.cpp"
#include "readmesh.cpp"
#include "makemesh.cpp"
#include "makemaster.cpp"
#include "domaindecomposition.cpp"
#include "connectivity.cpp"
#include "writebinaryfiles.cpp"
#include "CodeGenerator.cpp"
#include "CodeCompiler.cpp"

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        std::cerr << "Usage: ./parseinput <pdeapp.txt>\n";
        return 1;
    }    

    if (std::filesystem::exists(argv[1])) {
        std::cout << "Generating Exasim's input files for this text file ("<< argv[1] << ") ... \n\n";
    } else {
        error("Error: Input file does not exist.\n");        
    }          
           
    InputParams params = parseInputFile(argv[1]);                           
    PDE pde = initializePDE(params);    
        
    Mesh mesh = initializeMesh(params, pde);        
    Master master = initializeMaster(pde, mesh);                            
        
    writeBinaryFiles(pde, mesh, master);

#ifdef USE_CMAKE
    ParsedSpec spec = generateCppCode(pde);
#else    
    ParsedSpec spec = generateCppCode(pde);
    executeCppCode(spec); 
    buildDynamicLibraries(spec);     
    std::cout << "\n******** Done with generating input files and dynamic libraries for EXASIM ********\n";
#endif        
    
    return 0;
}
