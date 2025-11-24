#include "exasim.hpp"
#include "BuiltIn/Poisson3D/model.hpp"
#include "libexasim.h"
#include <algorithm>
#include <vector>

char** remove_first_arg(int argc, char** argv) {
    std::vector<char*> v;

    // Keep argv[0], skip argv[1], keep argv[2..]
    v.reserve(argc); 
    v.push_back(argv[0]);        // keep program name

    for (int i = 2; i < argc; ++i) {
        v.push_back(argv[i]);
    }

    v.push_back(nullptr);        // execv requires NULL termination

    // Allocate a real C array for execv-style calls
    char** new_argv = new char*[v.size()];
    std::copy(v.begin(), v.end(), new_argv);

    return new_argv;             // caller frees new_argv (NOT the strings)
}


enum Model {
  Poisson3DModel
};


int libmain(int argc, char ** argv){
  auto rest =  remove_first_arg(argc, argv);
  if (std::string(argv[1])  == "Poisson3D") {
  premain<Poisson3D::Poisson3D>(argc - 1, rest);
  }
  else {
    error("Error: Model does not exist!");
  }
}
