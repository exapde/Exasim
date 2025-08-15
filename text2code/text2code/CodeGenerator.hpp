/**
 * @class CodeGenerator
 * @brief Generates C++ code files based on parsed specifications.
 *
 * This class provides methods to generate various C++ source and header files
 * from a given ParsedSpec object. It supports code generation for symbolic functions,
 * scalars, vectors, CUDA/HIP headers, and several empty template files for different
 * model components. The generated files are intended to facilitate the implementation
 * of PDE models and related computational routines.
 *
 * @note All file generation methods take a filename or model path as input and write
 *       the corresponding code to the specified location.
 *
 * @see ParsedSpec
 * @see FunctionDef
 */

#pragma once
#include "TextParser.hpp"

class CodeGenerator {
public:
    CodeGenerator(const ParsedSpec& spec);

    void generateCode2Cpp(const std::string& filename) const;
    void generateSymbolicFunctionsHpp(const std::string& filename) const;
    void generateSymbolicFunctionsCpp(const std::string& filename) const;
    void generateSymbolicScalarsVectorsHpp(const std::string& filename) const;
    void generateSymbolicScalarsVectorsCpp(const std::string& filename) const;
    void generateCudaHipHpp(const std::string& filename) const;
        
//     void generateEmptyFluxCpp(std::string modelpath) const;
//     void generateEmptySourceCpp(std::string modelpath) const;
//     void generateEmptyTdfuncCpp(std::string modelpath) const;
//     void generateEmptyUbouCpp(std::string modelpath) const;
//     void generateEmptyFbouCpp(std::string modelpath) const;
//     void generateEmptyFbouHdgCpp(std::string modelpath) const;
    void generateEmptySourcewCpp(std::string modelpath) const;
    void generateEmptyOutputCpp(std::string modelpath) const;
    void generateEmptyMonitorCpp(std::string modelpath) const;
    void generateEmptyInituCpp(std::string modelpath) const;
    void generateEmptyInitqCpp(std::string modelpath) const;
    void generateEmptyInitudgCpp(std::string modelpath) const;
    void generateEmptyInitwdgCpp(std::string modelpath) const;
    void generateEmptyInitodgCpp(std::string modelpath) const;
    void generateEmptyAvfieldCpp(std::string modelpath) const;
    void generateEmptyFintCpp(std::string modelpath) const;    
    void generateEmptyFhatCpp(std::string modelpath) const;    
    void generateEmptyUhatCpp(std::string modelpath) const;    
    void generateEmptyStabCpp(std::string modelpath) const;    
    void generateEmptyEoSCpp(std::string modelpath) const;    
    void generateLibPDEModelHpp(std::string modelpath) const;
    void generateLibPDEModelCpp(std::string modelpath) const;
    
private:
    const ParsedSpec& spec;

    void generateFunctionHeader(std::ostream& os, const FunctionDef& func) const;
    void generateFunctionSource(std::ostream& os, const FunctionDef& func) const;
};




