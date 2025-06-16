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
        
//     void generateEmptyFluxCpp() const;
//     void generateEmptySourceCpp() const;
//     void generateEmptyTdfuncCpp() const;
//     void generateEmptyUbouCpp() const;
//     void generateEmptyFbouCpp() const;
//     void generateEmptyFbouHdgCpp() const;
    void generateEmptySourcewCpp() const;
    void generateEmptyOutputCpp() const;
    void generateEmptyMonitorCpp() const;
    void generateEmptyInituCpp() const;
    void generateEmptyInitqCpp() const;
    void generateEmptyInitudgCpp() const;
    void generateEmptyInitwdgCpp() const;
    void generateEmptyInitodgCpp() const;
    void generateEmptyAvfieldCpp() const;
    void generateEmptyFintCpp() const;    
    void generateEmptyFhatCpp() const;    
    void generateEmptyUhatCpp() const;    
    void generateEmptyStabCpp() const;    
    void generateEmptyEoSCpp() const;    
    void generateLibPDEModelHpp() const;
    void generateLibPDEModelCpp() const;
    
private:
    const ParsedSpec& spec;

    void generateFunctionHeader(std::ostream& os, const FunctionDef& func) const;
    void generateFunctionSource(std::ostream& os, const FunctionDef& func) const;
};




