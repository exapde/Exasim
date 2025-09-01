/*
 * text2code.cpp
 *
 * This program generates input files and dynamic libraries for the EXASIM framework
 * from a user-provided PDE application specification in a text file.
 *
 * Usage:
 *   ./parseinput <pdeapp.txt>
 *
 * Main Steps:
 *   1. Checks if the input file exists.
 *   2. Parses the input file to extract PDE and mesh parameters.
 *   3. Initializes PDE, mesh, and master structures.
 *   4. Writes binary files required by EXASIM.
 *   5. Generates C++ code for the PDE specification.
 *   6. Optionally compiles and builds dynamic libraries for EXASIM.
 *
 * Dependencies:
 *   - C++17 standard library
 *   - BLAS and LAPACK libraries
 *   - METIS and GKlib (optional, if HAVE_METIS is defined)
 *   - Several local source files: tinyexpr.cpp, helpers.cpp, readpdeapp.cpp, readmesh.cpp,
 *     makemesh.cpp, makemaster.cpp, domaindecomposition.cpp, connectivity.cpp,
 *     writebinaryfiles.cpp, CodeGenerator.cpp, CodeCompiler.cpp
 *
 * Compilation Examples:
 *   See commented lines at the top of the file for various compilation options.
 *
 * Author: Ngoc Cuong Nguyen
 * Date: 07/07/2025
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <array>
#include <regex>
#include <numeric>
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

#include "TextParser.hpp"

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
     
    ParsedSpec spec = TextParser::parseFile(make_path(pde.datapath, pde.modelfile));        
    spec.exasimpath = pde.exasimpath;
    spec.modelpath = make_path(pde.exasimpath, "/backend/Model/");
    spec.symenginepath = make_path(pde.exasimpath, "/text2code/symengine");
    
    Mesh mesh = initializeMesh(params, pde);        
    Master master = initializeMaster(pde, mesh);                            
        
    writeBinaryFiles(pde, mesh, master, spec);
    
#ifdef USE_CMAKE
    if (pde.gencode==1) generateCppCode(spec);
#else    
    if (pde.gencode==1) {
        generateCppCode(spec);
        executeCppCode(spec); 
        buildDynamicLibraries(spec);     
    }
#endif        
    
    std::cout << "\n******** Done with generating input files and dynamic libraries for EXASIM ********\n";

    return 0;
}
