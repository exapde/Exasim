#include <cstdlib>
#include <iostream>
#include <sstream>
#include <filesystem>
#include "CodeGenerator.cpp"

// g++ -O2 -std=c++17 Text2Code.cpp -o text2code

namespace fs = std::filesystem;

void moveSourceFiles(const fs::path& sourceDir, const fs::path& targetDir) {
    if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir)) {
        std::cerr << "Source directory does not exist.\n";
        return;
    }

    if (!fs::exists(targetDir)) {
        fs::create_directories(targetDir);
    }

    for (const auto& entry : fs::directory_iterator(sourceDir)) {
        if (entry.is_regular_file()) {
            fs::path file = entry.path();
            if (file.extension() == ".cpp" || file.extension() == ".h" || file.extension() == ".hpp") {
                fs::path targetPath = targetDir / file.filename();
                fs::rename(file, targetPath);  // Move the file
                std::cout << "Moved: " << file << " -> " << targetPath << "\n";
            }
        }
    }
}

std::string trimToSubstring(const std::string& fullPath, const std::string& keyword) {
    std::size_t pos = fullPath.find(keyword);  // Use rfind to get the last occurrence
    if (pos != std::string::npos) {
        return fullPath.substr(0, pos + keyword.length());
    }
    else {      
      return "";
    }
    return fullPath; // Return unmodified if keyword not found
}

// std::string getSharedLibExtension() {
// #if defined(_WIN32)
//     return ".dll";
// #elif defined(__APPLE__)
//     return ".dylib";
// #else // Linux and other POSIX systems
//     return ".so";
// #endif
// }

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./text2code <pdemodel.txt> [<text2code_path>]\n";
        return 1;
    }    
               
    // get current working directory
    std::filesystem::path cwd = std::filesystem::current_path();
            
    ParsedSpec spec;    
    
    try {
        if (std::filesystem::exists(argv[1])) {
            std::cout << "Generating the C++ code for this input text file ("<< argv[1] << ") ... \n\n";
        } else {
            std::cout << "Input file does not exist.\n";
            return 1;
        }  
        
        // parse the input text file 
        spec = TextParser::parseFile(argv[1]);

        // Generate SymEngine code
        CodeGenerator gen(spec);
        gen.generateCode2Cpp("Code2Cpp.cpp");
        gen.generateSymbolicFunctionsHpp("SymbolicFunctions.hpp");
        gen.generateSymbolicFunctionsCpp("SymbolicFunctions.cpp");        
        gen.generateSymbolicScalarsVectorsHpp("SymbolicScalarsVectors.hpp");
        gen.generateSymbolicScalarsVectorsCpp("SymbolicScalarsVectors.cpp");   
        
        std::string framework = spec.framework;
        if ((framework == "Cuda") || (framework == "CUDA") || (framework == "cuda")) {
          gen.generateCudaHipHpp("CudaHip.hpp");
        } 
        else if ((framework == "Hip") || (framework == "HIP") || (framework == "hip")) {
          gen.generateCudaHipHpp("CudaHip.hpp");
        } 
        
        if (spec.exasim == true) {
          if (spec.isoutput[6]==false) {gen.generateEmptySourcewCpp();}
          if (spec.isoutput[7]==false) {gen.generateEmptyOutputCpp();}
          if (spec.isoutput[8]==false) {gen.generateEmptyMonitorCpp();}
          if (spec.isoutput[9]==false) {gen.generateEmptyInituCpp();}
          if (spec.isoutput[10]==false) {gen.generateEmptyInitqCpp();}
          if (spec.isoutput[11]==false) {gen.generateEmptyInitudgCpp();}
          if (spec.isoutput[12]==false) {gen.generateEmptyInitwdgCpp();}
          if (spec.isoutput[13]==false) {gen.generateEmptyInitodgCpp();}
          if (spec.isoutput[14]==false) {gen.generateEmptyAvfieldCpp();}
          if (spec.isoutput[15]==false) {gen.generateEmptyFintCpp();}
          if (spec.isoutput[16]==false) {gen.generateEmptyEoSCpp();}
          gen.generateEmptyFhatCpp();
          gen.generateEmptyUhatCpp();
          gen.generateEmptyStabCpp();
          gen.generateLibPDEModelHpp();
          gen.generateLibPDEModelCpp();
        }
        
        std::cout << "The C++ code is generated in the following files:\n";
        std::cout << "    Code2Cpp.cpp\n";
        std::cout << "    SymbolicFunctions.hpp and SymbolicFunctions.cpp\n";
        std::cout << "    SymbolicScalarsVectors.hpp and SymbolicScalarsVectors.cpp\n";
               
        std::cout << "These files are located in this directory: " << cwd << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error during parsing or generation: " << ex.what() << "\n";
        return 1;
    }
    
    // get path to text2code   
    std::string text2code_path;
    if (argc >= 3) {
      text2code_path = argv[2];
    }
    else {
      text2code_path = trimToSubstring(cwd, "text2code");
      if (text2code_path == "") {
        text2code_path = trimToSubstring(cwd, "Exasim");            
        if (text2code_path == "") {
          std::cout <<"The text2code path could not be found. Please input the correct path to text2code when you run text2code.\n";
          return 1;
        }
        else { text2code_path = text2code_path + "/text2code"; }
      }      
    }
        
    if (std::filesystem::exists(text2code_path) && std::filesystem::is_directory(text2code_path)) {
        std::cout << "\nCompiling the C++ generated code (Code2Cpp.cpp) ...\n";
    } else {
        std::cout <<"This text2code path ("<<text2code_path<<") is not valid. Please open Text2Code.cpp and replace text2code_path with the correct path to text2code.\n";
        return 1;
    }
                
    // Construct compile command using text2code_path
    std::stringstream cmd;
    cmd << "g++ -std=c++17 -Wno-inconsistent-missing-override Code2Cpp.cpp "        
        //<< "-include " << text2code_path << "/text2code/pch.hpp "           
        << "-I" << text2code_path << " "
        << "-I" << text2code_path << "/include "       
        << text2code_path << "/lib/libsymengine.a "    
        << "-o code2cpp";    
        //<< "-L" << text2code_path << "/lib "        
        //<< "-Wl,-rpath," << text2code_path << "/lib "        
        //<< "-lsymengine -o code2cpp";
                   
    // g++ -std=c++17  -Wno-inconsistent-missing-override code2cpp.cpp -I/Users/ngoccuongnguyen/GitHub/text2code/ -I/Users/ngoccuongnguyen/GitHub/text2code/include -L/Users/ngoccuongnguyen/GitHub/text2code/lib -lsymengine -o code2cpp    
    std::cout<<cmd.str().c_str()<<std::endl;
    
    int status = std::system(cmd.str().c_str());
    if (status != 0) {
        std::cerr << "Compilation failed!\n";
        return 1;
    }    
                
    std::cout << "\nRunning the C++ generated code (Code2Cpp) ...\n\n";
    status = std::system("./code2cpp");
    if (status != 0) {
        std::cerr << "Generator execution failed!\n";
        return 1;
    }
    else if (spec.exasim == false) {
        std::cout << "The below source codes are generated and stored in this directory: " << cwd << "\n";      
        for (const auto& funcname : spec.outputs) {
          std::cout<<funcname<<".h, "<<funcname<<".cpp\n"; 
        }
    }
    else if (spec.exasim == true) {
      
      std::string keyword = "Exasim";
      std::string exasim_path = cwd;
      std::size_t pos = exasim_path.rfind(keyword);  // Use rfind to get the last occurrence
      
      if (pos != std::string::npos) {
          exasim_path = exasim_path.substr(0, pos + keyword.length());
          std::cout << "\nThe below source codes are generated and moved to this directory: " << exasim_path + "/backend/Model" << "\n";      
          moveSourceFiles(cwd, exasim_path + "/backend/Model");
          
//           std::cout << "\nCompiling the above source codes to generate dynamic libraries  ...\n\n";
//           
//           std::string kokkos_path = exasim_path + "/kokkos";          
//           std::string kokkos_serial_path = exasim_path + "/kokkos/buildserial";
//           std::string kokkos_cuda_path = exasim_path + "/kokkos/buildcuda";
//           std::string kokkos_hip_path = exasim_path + "/kokkos/buildhip";
//           
//           if (std::filesystem::exists(kokkos_serial_path) && std::filesystem::is_directory(kokkos_serial_path)) {
//               std::stringstream cmp;
//               cmp << "g++ -fPIC -shared -std=c++17 "        
//                   << exasim_path << "/backend/Model/libpdemodel.cpp "   
//                   << "-I" << kokkos_serial_path << "/include "       
//                   << kokkos_serial_path << "/lib/libkokkoscore.a "    
//                   << kokkos_serial_path << "/lib/libkokkoscontainers.a "        
//                   << kokkos_serial_path << "/lib/libkokkossimd.a -o "        
//                   << exasim_path << "/backend/Model/libpdemodelserial"<< getSharedLibExtension();                                   
// 
//               std::cout<<cmp.str().c_str()<<std::endl;
//               status = std::system(cmp.str().c_str());
//               if (status != 0) {
//                 std::cerr << "Compilation failed!" << std::endl;
//                 return 1;
//               }            
//           }                   
//           
//           if (std::filesystem::exists(kokkos_cuda_path) && std::filesystem::is_directory(kokkos_cuda_path)) {
//               std::stringstream cmp;
//               cmp << "clang++ -fPIC -shared -std=c++17 "     
//                   << exasim_path << "/backend/Model/libpdemodel.cpp "       
//                   << "-I" << kokkos_cuda_path << "/include "       
//                   << kokkos_cuda_path << "/lib/libkokkoscore.a "    
//                   << kokkos_cuda_path << "/lib/libkokkoscontainers.a -o "                             
//                   << exasim_path << "/backend/Model/libpdemodelcuda"<< getSharedLibExtension();                      
// 
//               std::cout<<cmp.str().c_str()<<std::endl;
//               status = std::system(cmp.str().c_str());
//               if (status != 0) {
//                 std::cerr << "Compilation failed!" << std::endl;
//                 return 1;
//               }            
//           }                   
//           
//           if (std::filesystem::exists(kokkos_hip_path) && std::filesystem::is_directory(kokkos_hip_path)) {
//               std::stringstream cmp;
//               cmp << "clang++ -fPIC -shared -std=c++17 libpdemodel.cpp "        
//                   << "-I" << kokkos_hip_path << "/include "       
//                   << kokkos_hip_path << "/lib/libkokkoscore.a "    
//                   << kokkos_hip_path << "/lib/libkokkoscontainers.a -o "                             
//                   << exasim_path << "/backend/Model/libpdemodelhip"<< getSharedLibExtension();        
// 
//               std::cout<<cmp.str().c_str()<<std::endl;
//               status = std::system(cmp.str().c_str());
//               if (status != 0) {
//                 std::cerr << "Compilation failed!" << std::endl;
//                 return 1;
//               }            
//           }                             
      }
      else {
        std::cout <<"The source codes are generated and stored in this directory: " << exasim_path << "\n";      
        std::cout <<"The current working directory \""<<exasim_path<<"\" does not contain \""<<keyword<<"\"\n";
        std::cout <<"Please make sure you run \"text2code\" from a directory containing \""<<keyword<<"\"\n";
        std::cout <<"Please move the generated source codes from " << exasim_path << "to /Exasim/backend/Model \n";      
        return 1;
      }
          
    }        

    return 0;
}

