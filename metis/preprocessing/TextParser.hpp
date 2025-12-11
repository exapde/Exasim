#pragma once
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <regex>

// struct FunctionDef {
//     std::string name;
//     std::string output;
//     int outputsize = 0;
//     std::vector<std::string> args;
//     std::vector<std::string> body;    
//     std::unordered_map<std::string, std::pair<int, int>> matrices; // name -> (rows, cols)
// };
// 
// struct ParsedSpec {
//     std::vector<std::string> scalars;
//     std::unordered_map<std::string, int> vectors;
//     std::vector<std::string> namevectors;
//     std::vector<std::string> jacobian;
//     std::vector<std::string> hessian;
//     std::vector<std::string> batch;
//     std::vector<std::string> outputs;
//     std::vector<std::string> exasimfunctions = {
//         "Flux", "Source", "Tdfunc", "Ubou", "Fbou", "FbouHdg",
//         "Sourcew", "Output", "Monitor", "Initu", "Initq", "Inituq",
//         "Initw", "Initv", "Avfield", "Fint", "EoS", "VisScalars", 
//         "VisVectors", "VisTensors", "QoIvolume", "QoIboundary"};
//     std::vector<bool> isoutput;     
//     std::string datatype = "dstype";
//     std::string framework = "kokkos";
//     std::string codeformat = "exasim";    
//     std::string exasimpath = "";
//     std::string modelpath = "";
//     std::string symenginepath = "";
//     std::string modelfile = "";
//     bool exasim;
//     std::vector<FunctionDef> functions;
// };

class TextParser {
public:
    static ParsedSpec parseFile(const std::string& filename) {
      
        ParsedSpec spec;
        spec.modelfile = filename;
        
        std::ifstream infile(filename);
        std::string line;

        std::regex kv_re(R"((\w+)\s+(.+))");
        std::regex vector_re(R"((\w+)\((\d+)\))");
        std::regex func_re(R"(^\s*function\s+(\w+)\(([^)]*)\)\s*$)");
        std::regex matrix_decl_re(R"(^\s*matrix\s+(\w+)\((\d+),(\d+)\);?\s*$)");
        std::regex output_size_re(R"(^\s*output_size\((\w+)\)\s*=\s*(\d+)\s*;?\s*$)");

        bool inFunction = false;
        FunctionDef currentFunc;

        while (std::getline(infile, line)) {
            line = trim(line);
            if (line.empty()) continue;

            if (line.rfind("function", 0) == 0) {
                std::smatch match;
                if (std::regex_match(line, match, func_re)) {
                    if (inFunction) {
                        spec.functions.push_back(currentFunc);
                        currentFunc = FunctionDef();
                    }
                    inFunction = true;                   
                    currentFunc.name = match[1];
                    currentFunc.args = split(match[2], ',');
                }
            } else if (inFunction) {
                std::smatch match;                
                if (std::regex_match(line, match, matrix_decl_re)) {
                    std::string name = match[1];
                    int rows = std::stoi(match[2]);
                    int cols = std::stoi(match[3]);                    
                    currentFunc.matrices[name] = std::make_pair(rows, cols);
                    currentFunc.body.push_back(line);
                } else if (std::regex_match(line, match, output_size_re)) {
                    std::string vec_name = match[1];                    
                    int size = std::stoi(match[2]);
                    currentFunc.output = match[1];
                    currentFunc.outputsize = size;                    
                    currentFunc.body.push_back(line);
                } else {
                    currentFunc.body.push_back(line);
                }
            } else {
                std::smatch match;
                if (std::regex_match(line, match, kv_re)) {
                    std::string key = match[1];
                    std::string value = match[2];
                    auto tokens = split(value, ',');

                    if (key == "scalars") {
                        spec.scalars = tokens;                        
                    } else if (key == "vectors") {
                        for (const auto& tok : tokens) {
                            std::smatch m;
                            if (std::regex_match(tok, m, vector_re)) {
                                spec.vectors[m[1]] = std::stoi(m[2]);                                
                                spec.namevectors.push_back(m[1]);                                
                            }
                        }
                    } else if (key == "jacobian") {
                        spec.jacobian = tokens;
                    } else if (key == "hessian") {
                        spec.hessian = tokens;
                    } else if (key == "batch") {
                        spec.batch = tokens;
                    } else if (key == "outputs") {
                        spec.outputs = tokens;
                    } else if (key == "datatype") {
                        spec.datatype = value;
                    } else if (key == "framework") {
                        spec.framework = value;                        
                    } else if (key == "codeformat") {
                        spec.codeformat = value;    
                    }
                }
            }
        }

        if (inFunction) {
            spec.functions.push_back(currentFunc);
        }
        
//         for (int i=0; i<spec.functions.size(); i++)
//           if (spec.functions[i].outputsize == 0) {
//             printf("Error: Keyword output_size is missing in a function body\n");
//             exit(-1);            
//           }
                          
        if ((spec.codeformat == "exasim") || (spec.codeformat == "Exasim") || (spec.codeformat == "EXASIM")) {
          spec.exasim = true;
          std::unordered_set<std::string> outputSet(spec.outputs.begin(), spec.outputs.end());

          spec.isoutput.reserve(spec.exasimfunctions.size());          
          for (const auto& name : spec.exasimfunctions) {
              spec.isoutput.push_back(outputSet.count(name) > 0);
          }

          // spec.isoutput.resize(spec.exasimfunctions.size()); 
          // for (int i=0; i < spec.exasimfunctions.size(); i++) {
          //     std::string name = spec.exasimfunctions[i];
          //     spec.isoutput[i] = (outputSet.count(name) > 0) ? true : false;
          // }

          for (size_t i = 0; i < 6; ++i) {
            if (spec.isoutput[i] == false) {
              std::cerr << "Error: \"" << spec.exasimfunctions[i]
                        << "\" is not defined in the model file. Please define it.\n";
              std::exit(EXIT_FAILURE);  // More conventional than -1
            }
          }                    
        } else {
          spec.exasim = false;                        
        }
        
        // for (const auto& vec : spec.vectors) {
        //     const std::string& name = vec.first;
        //     int size = vec.second;
        //     std::cout<<name<<" : "<<size<<std::endl;
        // }
        
        return spec;
    }
    
private:
    static std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t");
        size_t end = s.find_last_not_of(" \t");
        return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    }

    static std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delimiter)) {
            tokens.push_back(trim(item));
        }
        return tokens;
    }
};
