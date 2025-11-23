/*
 * CodeGenerator.cpp
 * 
 * This file implements the CodeGenerator class, which is responsible for generating C++ source and header files
 * for symbolic PDE models using SymEngine. The generator supports multiple frameworks (C++, CUDA, HIP, Kokkos)
 * and produces code for model functions, their Jacobians, Hessians, and initialization routines.
 * 
 * Key functionalities:
 * - Symbolic function code generation (headers and sources)
 * - Automatic differentiation (Jacobian and Hessian code generation)
 * - Common Subexpression Elimination (CSE) for efficient code
 * - Support for vector and scalar symbolic inputs
 * - Generation of empty stubs for model functions (for extensibility)
 * - Framework-specific parallel loop emission (C++, CUDA, HIP, Kokkos)
 * - Generation of libpdemodel.hpp/cpp for model interface aggregation
 * 
 * Main methods:
 * - generateCode2Cpp: Generates main driver code for symbolic function evaluation.
 * - generateCudaHipHpp: Emits GPU backend abstraction header for CUDA/HIP.
 * - generateSymbolicFunctionsHpp/Cpp: Generates headers and sources for symbolic functions.
 * - generateSymbolicScalarsVectorsHpp/Cpp: Generates headers and sources for symbolic input management.
 * - generateEmpty*Cpp: Emits empty stubs for various model functions.
 * - generateLibPDEModelHpp/Cpp: Aggregates all model function interfaces and sources.
 * 
 * Helper methods:
 * - emitSymbolicScalarsVectors, emitevaluateSymbolicFunctions, emitfunc2cse, emitfuncjac2cse, emitfuncjachess2cse:
 *   Emit code for symbolic input initialization, function evaluation, and CSE routines.
 * - emitfunc2cppfiles, emitfuncjac2cppfiles, emitfuncjachess2cppfiles, emitinitfunc2cppfiles:
 *   Emit code for function, Jacobian, Hessian, and initialization routines.
 * - appendUbouFbou, appendFbouHdg: Emit code for boundary condition function aggregation.
 * 
 * Usage:
 * 1. Construct CodeGenerator with a ParsedSpec describing the model.
 * 2. Call the appropriate generate* methods to emit code files for the model.
 * 
 * Dependencies:
 * - SymEngine (symbolic computation)
 * - Standard C++ libraries (fstream, iostream, regex, etc.)
 * 
 * Author: [Your Name]
 * Date: [Date]
 */
#include "CodeGenerator.hpp"
#include <iostream>

std::string prefixSymEngineFunctions(const std::string& expr) {
    static const std::regex func_regex(
    //R"(\b(sin|cos|tan|cot|asin|acos|atan|acot|sinh|cosh|tanh|coth|asinh|acosh|atanh|acoth|sech|csch|log|log10|exp|sqrt|cbrt|abs|ceil|floor|pow|max|min|sign|pi)\b)");
            R"(\b(pi)\b)");
    return std::regex_replace(expr, func_regex, "SymEngine::$1");
}

CodeGenerator::CodeGenerator(const ParsedSpec& spec_) : spec(spec_) {}

void CodeGenerator::generateCode2Cpp(const std::string& filename) const {
    std::ofstream os(filename);

    os << "#include \"SymbolicFunctions.cpp\"\n";
    os << "#include \"SymbolicScalarsVectors.cpp\"\n\n";
                
    os << "int main() \n";
    os << "{\n";
    os << "  SymbolicScalarsVectors ssv;\n\n";
    os << "  for (int i=0; i<ssv.outputfunctions.size(); i++) {\n";
    os << "    std::string funcname = ssv.funcnames[i];\n";    
    os << "    if (ssv.outputfunctions[i] == true) {\n";
    os << "      std::vector<Expression> f = ssv.evaluateSymbolicFunctions(i);\n";
    os << "      bool append = true;\n";
    os << "      if (ssv.exasim == true) {\n";
    os << "        std::string fname = std::string(\"Kokkos\")  + funcname;\n";
    os << "        std::string jname = std::string(\"Hdg\") + funcname;\n";
    os << "        if (funcname == \"FbouHdg\") {\n"; 
    os << "          fname = std::string(\"HdgFbouonly\");\n";
    os << "          jname = std::string(\"HdgFbou\");\n";
    os << "        }\n";
    os << "        if (funcname == \"Sourcew\") {\n"; 
    os << "          fname = std::string(\"HdgSourcewonly\");\n";
    os << "          jname = std::string(\"HdgSourcew\");\n";
    os << "        }\n";    
    os << "        if (funcname == \"QoIboundary\") { \n";
    os << "            ssv.func2cppfiles(f, ssv.modelpath + fname, fname + std::to_string(1), i, false);\n";
    os << "            ssv.appendUbouFbou(ssv.modelpath + fname, fname, 1);\n";
    os << "        }\n";
    //os << "        std::cout<<funcname<<\", \"<<fname<<\", \"<<jname<<std::endl;\n";    
    os << "        else if ((funcname == \"Ubou\") || (funcname == \"Fbou\") || (funcname == \"FbouHdg\")) { \n";
    os << "          int szf = f.size();\n";
    os << "          int szuhat = ssv.szuhat;\n";
    os << "          int nbc = szf/szuhat;\n";
    //os << "        std::cout<<szf<<\", \"<<szuhat<<\", \"<<nbc<<std::endl;\n";    
    os << "          for (int n = 0; n < nbc; ++n) {\n";
    os << "            std::vector<Expression> g(szuhat);\n";        
    os << "            for (int m = 0; m < szuhat; ++m) {\n";
    os << "              g[m] = f[m + n * szuhat];\n";
    os << "            }\n";
    os << "            if (n==0) {\n";
    os << "               ssv.func2cppfiles(g, ssv.modelpath + fname, fname + std::to_string(n+1), i, false);\n";
    os << "               if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(g, ssv.modelpath + jname, jname + std::to_string(n+1), i, false);\n";
    os << "            } else {\n";
    os << "               ssv.func2cppfiles(g, ssv.modelpath + fname, fname + std::to_string(n+1), i, append);\n";
    os << "               if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(g, ssv.modelpath + jname, jname + std::to_string(n+1), i, append);\n";
    os << "            }\n";
    os << "            if (n==nbc-1) {\n";
    os << "               ssv.appendUbouFbou(ssv.modelpath + fname, fname, nbc);\n";
    os << "               if (funcname == \"FbouHdg\") ssv.appendFbouHdg(ssv.modelpath + jname, jname, nbc);\n";
    os << "            }\n";
    os << "          }\n";    
    os << "        } else if ((funcname == \"Initu\") || (funcname == \"Initq\") || (funcname == \"Inituq\") || (funcname == \"Initv\") || (funcname == \"Initw\")) { \n";
    os << "          std::string kname = funcname;\n";    
    os << "          if (funcname == \"Inituq\") kname = \"Initudg\";\n"; 
    os << "          if (funcname == \"Initv\") kname = \"Initvdg\";\n"; 
    os << "          if (funcname == \"Initw\") kname = \"Initwdg\";\n"; 
    os << "          ssv.initfunc2cppfiles(f, ssv.modelpath + \"cpu\" + kname, \"cpu\" + kname, i, false, 0);\n";
    os << "          ssv.initfunc2cppfiles(f, ssv.modelpath + \"Kokkos\" + kname, \"Kokkos\" + kname, i, false, 1);\n";
    os << "        } else {\n";
    os << "          ssv.func2cppfiles(f, ssv.modelpath + fname, fname, i, false);\n";
    os << "          if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, ssv.modelpath + jname, jname, i, false);\n";
    os << "        }\n";
    os << "      } else {\n";
    os << "        ssv.func2cppfiles(f, funcname, funcname, i, false);\n";
    os << "        if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, funcname, funcname, i, append);\n";
    os << "        if (ssv.hessianInputs[i].size() > 0) ssv.funcjachess2cppfiles(f, funcname, funcname, i, append);\n";
    os << "      }\n";
    os << "    }\n";
    os << "  }\n";
    os << "}\n";
    
    os.close();  
}

void CodeGenerator::generateCudaHipHpp(const std::string& filename) const {  
  std::ofstream os(filename);
  
  os << "#pragma once\n\n";
  os << "#include <string>\n\n";
  os << "#if defined(__HIPCC__)\n";
  os << "  #include <hip/hip_runtime.h>\n";
  os << "  #define GPU_DEVICE __device__\n";
  os << "  #define GPU_GLOBAL __global__\n";
  os << "  #define gpuDeviceSynchronize hipDeviceSynchronize\n";
  os << "  #define GPU_BACKEND \"HIP\"\n";
  os << "#elif defined(__CUDACC__)\n";
  os << "  #include <cuda_runtime.h>\n";
  os << "  #define GPU_DEVICE __device__\n";
  os << "  #define GPU_GLOBAL __global__\n";
  os << "  #define gpuDeviceSynchronize cudaDeviceSynchronize\n";
  os << "  #define GPU_BACKEND \"CUDA\"\n";
  os << "#else\n";
  os << "  #error \"GPU backend not recognized. Use nvcc or hipcc.\"\n";
  os << "#endif\n\n";
  os << "#define GPU_LAMBDA [=] GPU_DEVICE\n\n";
  os << "template <typename Functor>\n";
  os << "GPU_GLOBAL void lambda_kernel(Functor f, size_t N) {\n";
  os << "    size_t i = blockIdx.x * blockDim.x + threadIdx.x;\n";
  os << "    if (i < N) f(i);\n";
  os << "}\n\n";
  os << "template <typename Functor>\n";
  os << "void parallel_for(const std::string& label, size_t N, Functor f, int blockSize = 256) {\n";
  os << "    int numBlocks = (N + blockSize - 1) / blockSize;\n";
  os << "    lambda_kernel<<<numBlocks, blockSize>>>(f, N);\n";
  os << "    gpuDeviceSynchronize();\n";
  os << "}\n";
  
  os.close();  
}

void CodeGenerator::generateSymbolicFunctionsHpp(const std::string& filename) const {
    std::ofstream os(filename);
    os << "#pragma once\n\n";
    
    os << "#include <vector>\n";
    os << "#include <string>\n";
    os << "#include <fstream>\n";
    os << "#include <iostream>\n";
    os << "#include <symengine/expression.h>\n";
    os << "// #include <symengine/functions.h>\n";
    os << "#include \"SymEngineFunctionWrappers.hpp\"\n";
    os << "#include <symengine/printers/codegen.h>\n";
    //os << "#include <symengine/cse.h>\n";
    os << "#include <symengine/symbol.h>\n";
    os << "#include <symengine/matrix.h>\n\n";

    os << "using SymEngine::Expression;\n";
    os << "using SymEngine::vec_pair;\n";
    os << "using SymEngine::vec_basic;\n";
    os << "using SymEngine::symbol;\n";
    os << "using SymEngine::RCP;\n";
    os << "using SymEngine::Basic;\n";
    os << "using SymEngine::map_basic_basic;\n";
    os << "using SymEngine::CodePrinter;\n";
    os << "using SymEngine::C99CodePrinter;\n";
    
    os << "\n";
    for (const auto& func : spec.functions) {
        generateFunctionHeader(os, func);
    }
    
    os.close();  
}

void CodeGenerator::generateSymbolicFunctionsCpp(const std::string& filename) const {
    std::ofstream os(filename);
    os << "#include \"SymbolicFunctions.hpp\"\n\n";
    
    for (const auto& func : spec.functions) {
        generateFunctionSource(os, func);
    }
    
    os.close();  
}

void CodeGenerator::generateFunctionHeader(std::ostream& os, const FunctionDef& func) const {
    os << "std::vector<Expression> " << func.name << "(";
    for (size_t i = 0; i < func.args.size(); ++i) {
        bool isscalar = false;       
        for (int j=0; j<spec.scalars.size(); j++) {               
          if (spec.scalars[j] == func.args[i]) {
            isscalar = true;
            break;
          }
        }
        bool isvector = false;       
        for (int j=0; j<spec.namevectors.size(); j++) {               
          if (spec.namevectors[j] == func.args[i]) {
            isvector = true;
            break;
          }
        }
        
        if ((isscalar==false) && (isvector==false)) {
          printf("Error: function argument %s is not listed in scalars or vectors\n", func.args[i].c_str());
          exit(-1);
        }
        
        if (isscalar==true) os << "const Expression& " << func.args[i];
        if (isvector==true) os << "const std::vector<Expression>& " << func.args[i];           
      
        if (i + 1 < func.args.size()) os << ", ";
    }
    os << ");\n";
}

void CodeGenerator::generateFunctionSource(std::ostream& os, const FunctionDef& func) const {
    std::string lhs = func.output;
    int lhs_size = func.outputsize;
    os << "std::vector<Expression> " << func.name << "(";
    for (size_t i = 0; i < func.args.size(); ++i) {
        bool isscalar = false;       
        for (int j=0; j<spec.scalars.size(); j++) {               
          if (spec.scalars[j] == func.args[i]) {
            isscalar = true;
            break;
          }
        }
        bool isvector = false;       
        for (int j=0; j<spec.namevectors.size(); j++) {               
          if (spec.namevectors[j] == func.args[i]) {
            isvector = true;
            break;
          }
        }
        
        if ((isscalar==false) && (isvector==false)) {
          printf("Error: function argument %s is not listed in scalars or vectors\n", func.args[i].c_str());
          exit(-1);
        }
        
        if (isscalar==true) os << "const Expression& " << func.args[i];
        if (isvector==true) os << "const std::vector<Expression>& " << func.args[i];           
        
        if (i + 1 < func.args.size()) os << ", ";
    }
    os << ") {\n";

    os << "    std::vector<Expression> " << lhs << ";\n";
    os << "    " << lhs << ".resize(" << lhs_size << ");\n\n";
    
    std::regex for_loop_start(R"(^\s*for\s+(\w+)\s+in\s+(\d+):(\d+)\s*$)");
    std::regex for_loop_end(R"(^\s*endfor\s*$)");
    bool in_loop = false;

    std::regex call_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\((.*)\)\s*$)");
    std::regex assign_pattern(R"(^\s*(\w+)\s*=\s*(.+)$)");    
    std::regex zeros_pattern(R"(^\s*zeros\s*\(\s*(\w+)\s*,\s*(\d+)\s*\)\s*$)");
    std::regex  ones_pattern(R"(^\s*ones\s*\(\s*(\w+)\s*,\s*(\d+)\s*\)\s*$)");
    //std::regex fill_pattern(R"(^\s*fill\((\w+),\s*(\d+)\s*,\s*(\d+)\)\s*;?\s*$)");
    std::regex fill_pattern(R"(^\s*fill\((\w+),\s*(\d+),\s*([-+]?[0-9]*\.?[0-9]+)\)\s*;?\s*$)");
    std::regex zeros_no_size_pattern(R"(^\s*zeros\((\w+)\)\s*;?\s*$)");
    std::regex ones_no_size_pattern(R"(^\s*ones\((\w+)\)\s*;?\s*$)");
    std::regex fill_no_size_pattern(R"(^\s*fill\((\w+),\s*([-+]?[0-9]*\.?[0-9]+)\)\s*;?\s*$)");

    std::regex vec_decl_pattern(R"(^\s*vector\s+(\w+)\((\d+)\)\s*$)");
    
    std::regex mat_decl_pattern(R"(^\s*matrix\s+(\w+)\((\d+),(\d+)\)\s*$)");
    std::regex matset_pattern(R"((\w+)\[(\d+)\]\[(\d+)\]\s*=\s*(.+))");
    std::regex matset_pattern_ext(R"((\w+)\[([^\]]+)\]\[([^\]]+)\]\s*=\s*(.+))");
    std::regex matget_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\[(\d+)\]\[(\d+)\]\s*$)");
    std::regex matget_pattern_ext(R"(^\s*(.+)\s*=\s*(\w+)\[([^\]]+)\]\[([^\]]+)\]\s*$)");
    std::regex det_pattern(R"(^\s*(\w+)\[(\d+)\]\s*=\s*det\((\w+)\)\s*$)");
    std::regex trace_pattern(R"(^\s*(\w+)\[(\d+)\]\s*=\s*trace\((\w+)\)\s*$)");
    std::regex inv_pattern(R"(^\s*(\w+)\s*=\s*inv\((\w+)\)\s*$)");
    std::regex transpose_pattern(R"(^\s*(\w+)\s*=\s*transpose\((\w+)\)\s*$)");
    std::regex mat_binop_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*([\+\*])\s*(\w+)\s*$)");
    std::regex scalar_mat_mul_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*\*\s*(\w+)\s*$)");
    std::regex matmul_pattern(R"(^\s*(\w+)\s*=\s*(\w+)\s*\*\s*(\w+)(\s*\*\s*(\w+))?\s*$)");

    for (const auto& line_raw : func.body) {
        std::string line = prefixSymEngineFunctions(line_raw);             
        
        line.erase(std::remove(line.begin(), line.end(), ';'), line.end());

        std::smatch match;

        if (std::regex_match(line, match, for_loop_start)) {
            std::string var = match[1];
            std::string start = match[2];
            std::string end = match[3];
            os << "    for (int " << var << " = " << start << "; " << var << " <= " << end << "; ++" << var << ") {\n";
            in_loop = true;
        } else if (std::regex_match(line, match, for_loop_end)) {
            os << "    }\n";
            in_loop = false;
        } else if (std::regex_match(line, match, vec_decl_pattern)) {
            std::string vec = match[1];
            std::string size = match[2];
            os << "    std::vector<Expression> " << vec << "(" << size << ");\n";    
        } else if (std::regex_match(line, match, mat_decl_pattern)) {
            std::string mat = match[1];
            std::string rows = match[2];
            std::string cols = match[3];
            os << "    SymEngine::DenseMatrix " << mat << "(" << rows << ", " << cols << ");\n";
        } else if (std::regex_match(line, match, matset_pattern)) {
            std::string mat = match[1];
            std::string i = match[2];
            std::string j = match[3];
            std::string rhs = match[4];
            os << "    " << mat << ".set(" << i << ", " << j << ", Expression(" << rhs << "));\n";
        } else if (std::regex_match(line, match, matset_pattern_ext)) {
            std::string mat = match[1];
            std::string i = match[2];
            std::string j = match[3];
            std::string rhs = match[4];
            os << "    " << mat << ".set(" << i << ", " << j << ", Expression(" << rhs << "));\n";            
        } else if (std::regex_match(line, match, matget_pattern)) {
            std::string lhs = match[1];
            std::string mat = match[2];
            std::string i = match[3];
            std::string j = match[4];
            os << "    Expression " << lhs << " = " << mat << ".get(" << i << ", " << j << ");\n";    
        } else if (std::regex_match(line, match, matget_pattern_ext)) {
            std::string lhs = match[1];
            std::string mat = match[2];
            std::string i = match[3];
            std::string j = match[4];
            os << "    " << lhs << " = " << mat << ".get(" << i << ", " << j << ");\n";    
        } else if (std::regex_match(line, match, det_pattern)) {
            std::string vec = match[1];
            std::string index = match[2];
            std::string mat = match[3];
            os << "    " << vec << "[" << index << "] = " << mat << ".det();\n";
        } else if (std::regex_match(line, match, trace_pattern)) {
            std::string vec = match[1];
            std::string index = match[2];
            std::string mat = match[3];
            os << "    " << vec << "[" << index << "] = " << mat << ".trace();\n";
        } else if (std::regex_match(line, match, inv_pattern)) {
            std::string lhs_mat = match[1];
            std::string rhs_mat = match[2];
            //std::cout<<match[1]<<", "<<match[2]<<"\n";
            //std::cout<<func.matrices.count(lhs_mat)<<"\n";
            //for (const auto& [name, dims] : func.matrices) {
            //  std::cout << name << " : " << dims.first << " x " << dims.second << std::endl;
            //}
            if (func.matrices.count(lhs_mat) && func.matrices.count(rhs_mat)) {
                //os << "    " << lhs_mat << " = " << rhs_mat << ".inv();\n";
                os << "    " << rhs_mat << ".inv(" << lhs_mat << ");\n";
            }
        } else if (std::regex_match(line, match, transpose_pattern)) {
            std::string lhs_mat = match[1];
            std::string rhs_mat = match[2];
            if (func.matrices.count(lhs_mat) && func.matrices.count(rhs_mat)) {
                //os << "    " << lhs_mat << " = " << rhs_mat << ".transpose();\n";
                os << "    " << rhs_mat << ".transpose(" << lhs_mat << ");\n";
            }
        } else if (std::regex_match(line, match, mat_binop_pattern)) {
            std::string lhs_mat = match[1];
            std::string mat1 = match[2];
            std::string op = match[3];
            std::string mat2 = match[4];
            if (func.matrices.count(lhs_mat) && func.matrices.count(mat1) && func.matrices.count(mat2)) {
                //os << "    " << lhs_mat << " = DenseMatrix(" << mat1 << ") " << op << " DenseMatrix(" << mat2 << ");\n";
                //os << "    " << lhs_mat << " = " << mat1 << " " << op << " " << mat2 << ";\n";
              if (op == "*") os << "    " << mat1 << ".mul_matrix("<< mat2 << ", " << lhs_mat <<");\n"; 
              else if (op == "+") os << "    " << mat1 << ".add_matrix("<< mat2 << ", " << lhs_mat <<");\n"; 
            } else {
                os << "    Expression " << lhs_mat << " = " << mat1 << " " << op << " " << mat2 << ";\n";
            }
        } else if (std::regex_match(line, match, scalar_mat_mul_pattern)) {
            std::string lhs = match[1];
            std::string scalar = match[2];
            std::string mat = match[3];
            if (func.matrices.count(mat)) {
                os << "    " << lhs << " = DenseMatrix(" << mat << ") * " << scalar << ";\n";
            } else {
                os << "    Expression " << lhs << " = " << scalar << " * " << mat << ";\n";
            }
        } else if (std::regex_match(line, match, matmul_pattern)) {
            std::string lhs = match[1];
            std::string m1 = match[2];
            std::string m2 = match[3];
            std::string m3 = match[5];
            bool lhs_is_matrix = func.matrices.count(lhs);
            bool rhs_has_matrix = func.matrices.count(m1) || func.matrices.count(m2) || (!m3.empty() && func.matrices.count(m3));
            if (lhs_is_matrix || rhs_has_matrix) {
                os << "    " << lhs << " = " << m1 << " * " << m2;
                if (!m3.empty()) os << " * " << m3;
                os << ";\n";
            } else {
                os << "    Expression " << lhs << " = " << m1 << " * " << m2;
                if (!m3.empty()) os << " * " << m3;
                os << ";\n";
            }
        } else if (std::regex_match(line, match, call_pattern)) {
            std::string var = match[1];
            std::string fname = match[2];
            std::string args = match[3];
            if (var != lhs)
              os << "    auto " << var << " = " << fname << "(" << args << ");\n";
            else
              os << "    " << var << " = " << fname << "(" << args << ");\n";
        } else if (std::regex_match(line, match, zeros_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression(0));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(0);\n";              
              os << "    }\n";
            }
        } else if (std::regex_match(line, match, zeros_no_size_pattern)) {
            std::string var = match[1];
            if (var != lhs) {
              printf("Error: zeros(%s) is invalid", var.c_str());
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(0);\n";              
              os << "    }\n";
            }            
         } else if (std::regex_match(line, match, ones_no_size_pattern)) {
            std::string var = match[1];
            if (var != lhs) {
              printf("Error: ones(%s) is invalid", var.c_str());
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(1);\n";              
              os << "    }\n";
            }                
        } else if (std::regex_match(line, match, ones_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression(1));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression(1);\n";              
              os << "    }\n";
            }
        } else if (std::regex_match(line, match, fill_pattern)) {
            std::string var = match[1];
            int size = std::stoi(match[2]);
            double value = std::stod(match[3]);
            if (var != lhs) {
              os << "    std::vector<Expression> " << var << "(" << size << ", Expression("<< value << "));\n";              
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression("<< value <<");\n";              
              os << "    }\n";
            }            
        } else if (std::regex_match(line, match, fill_no_size_pattern)) {
            std::string var = match[1];
            double value = std::stod(match[2]);
            if (var != lhs) {
              printf("Error: fill(%s, %g) is invalid", var.c_str(), value);
              exit(-1);
            }
            else {
              std::string iv = "i";
              os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << lhs_size << "; ++" << iv << ") {\n";
              os << "         " << var << "[" << iv << "] = Expression("<< value <<");\n";              
              os << "    }\n";
            }
        } else if (line.find(lhs + "[") == 0) {
            size_t eq = line.find('=');
            std::string lhs_expr = line.substr(0, eq);
            std::string rhs_expr = line.substr(eq + 1);
            os << "    " << lhs_expr << " = " << rhs_expr << ";\n";
        } else if (std::regex_match(line, match, assign_pattern)) {
            std::string var = match[1];
            std::string rhs = match[2];
            os << "    Expression " << var << " = " << rhs << ";\n";
        }
    }

    os << "    return " << lhs << ";\n";
    os << "}\n\n";
}

void CodeGenerator::generateSymbolicScalarsVectorsHpp(const std::string& filename) const {
    std::ofstream os(filename);
    if (!os) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    os << "#pragma once\n\n";
    os << "#include \"SymbolicFunctions.hpp\"\n\n";
    
    os << "class SymbolicScalarsVectors {\n\n";
    os << "public:\n\n";
    
    os << "    // path to model folder \n";
    os << "    std::string modelpath = \"" << spec.modelpath << "\";\n\n";
    
    // Scalars
    os << "    // input symbolic scalars\n";
    for (const auto& s : spec.scalars) {
        os << "    Expression " << s << ";\n";
    }
    os << "\n";

    // Vectors
    os << "    // input symbolic vectors\n";
    for (const auto& [name, size] : spec.vectors) {
        os << "    std::vector<Expression> " << name << ";\n";
    }
    os << "\n";

    // Vector sizes
    os << "    // vector sizes\n";
    for (const auto& [name, size] : spec.vectors) {
        os << "    int sz" << name << ";\n";
    }
    os << "    bool exasim;\n";
    os << "\n";
        
    os << "    std::vector<bool> outputfunctions;\n";
    os << "    std::vector<std::vector<std::string>> funcargs;\n";
    os << "    std::vector<std::vector<std::string>> funcargssizes;\n";
    os << "    std::vector<std::string> funcnames;\n";
    os << "    std::vector<std::string> funcdecls;\n";
    os << "    std::vector<std::string> funcjacdecls;\n\n";
    
    os << "    std::vector<std::vector<std::pair<std::string, std::vector<Expression>>>> inputvectors;\n";  
    os << "    std::vector<std::vector<std::pair<std::string, Expression>>> inputscalars;\n";  
    os << "    std::vector<std::vector<std::vector<Expression>>> jacobianInputs;\n";
    os << "    std::vector<std::vector<std::vector<Expression>>> hessianInputs;\n\n";
    os << "    std::vector<std::string> batch;\n";
    
    os << "    SymbolicScalarsVectors();\n\n";
    
    os << "    std::vector<Expression> evaluateSymbolicFunctions(int call);\n\n";    
        
    os << "    void func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f);\n\n";
    
    os << "    void funcjac2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n";
    os << "                         std::vector<vec_basic> &reduced_exprs_J, const std::vector<Expression> &f,\n";
    os << "                         const std::vector<std::vector<Expression>>& inputs_J);\n\n";

    os << "    void funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n";
    os << "                         std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,\n";
    os << "                         const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,\n";
    os << "                         const std::vector<std::vector<Expression>>& inputs_H);\n\n";
    
    os << "    void func2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);\n";
    os << "    void funcjac2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);\n";
    os << "    void funcjachess2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);\n";
    
    os << "    void initfunc2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append, int framework);\n";
    os << "    void appendUbouFbou(const std::string& filename, const std::string& funcname, int nbc);\n";
    os << "    void appendFbouHdg(const std::string& filename, const std::string& funcname, int nbc);\n";
    
    os << "};\n";
    
    os.close();  
}

void emitSymbolicScalarsVectors(std::ostream& os, const ParsedSpec& spec) {
  
    os << "SymbolicScalarsVectors::SymbolicScalarsVectors() {\n\n";

    for (const auto& s : spec.scalars) {
        os <<"    "<<s<< " =   Expression(\"" << s << "\")" << ";\n";
    }
    os << "\n";
                
    for (const auto& vec : spec.vectors) {
        const std::string& name = vec.first;
        int size = vec.second;
        std::string szname = "sz" + name;

        os << "    " << szname << " = " << size << ";\n";
        os << "    " << name << ".resize(" << size << ");\n";
        
        std::string iv = "i";
        os << "    for (int " << iv << " = " << 0 << "; " << iv << " < " << size << "; ++" << iv << ") {\n";
        os << "         " << name << "[" << iv << "] = Expression(\"" << name <<"\"  + std::to_string(i));\n";
        os << "    }\n";
                      
        os << "\n";
    }
    if (spec.exasim == true) os << "    exasim = true;\n\n";
    else os << "    exasim = false;\n\n";

    os << "    outputfunctions.assign(" << spec.functions.size() << ", false);\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
        const std::string& funcname = spec.functions[i].name;
        if (std::find(spec.outputs.begin(), spec.outputs.end(), funcname) != spec.outputs.end()) {
            os << "    outputfunctions[" << i << "] = true;\n";
        }
    }
    os << "\n";

    // Emit funcnames
    os << "    batch = {";
    for (size_t i = 0; i < spec.batch.size(); ++i) {
        os<<"\""<<spec.batch[i]<<"\"";
        if (i < spec.batch.size() - 1) os << ", ";
    }
    os << "};\n\n";
    
    // Emit funcnames
    os << "    funcnames = {";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
        os<<"\""<<spec.functions[i].name<<"\"";
        if (i < spec.functions.size() - 1) os << ", ";
    }
    os << "};\n\n";
    
    // Emit funcargs
    os << "    funcargs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";
        for (size_t i = 0; i < func.args.size(); ++i) {
            os << "\"" << func.args[i] << "\"";
            if (i < func.args.size() - 1) os << ", ";
        }
        os << "}";
        if (fi < spec.functions.size() - 1) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    // Emit funcargs
    os << "    funcargssizes = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";
        for (size_t i = 0; i < func.args.size(); ++i) {
            os << "\"sz" << func.args[i] << "\"";
            if (i < func.args.size() - 1) os << ", ";
        }
        os << "}";
        if (fi < spec.functions.size() - 1) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    funcdecls = {\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
      const FunctionDef& func = spec.functions[i];
      const std::string& type = spec.datatype;
      const std::string& name = func.name;
            
      os <<"       \""<<"void " << name << "(" << type << "* f";
      // Function arguments
      for (const std::string& arg : func.args) {
          bool is_scalar = std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end();
          os << ", const " << type;
          if (!is_scalar) os << "*";
          os << " " << arg;
      }
      
      // Metadata arguments
      if (spec.exasim == true)
        os << ", const int modelnumber, const int N";
      else
        os << ", const int N, const int szf";
      for (const std::string& arg : func.args) {
          if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
              os << ", const int sz" << arg;
          }
      }
      os << ")"<<"\"";
      if (i < spec.functions.size() - 1) os << ", ";
      os << "\n";
    }
    os << "    };\n\n";
                
    os << "    funcjacdecls = {\n";
    for (size_t i = 0; i < spec.functions.size(); ++i) {
      const FunctionDef& func = spec.functions[i];
      const std::string& type = spec.datatype;
      const std::string& name = func.name;
            
      os <<"       \"";
      // Function arguments
      for (int i=0; i < func.args.size(); i++) {
          std::string arg = func.args[i];
          bool is_scalar = std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end();
          if (i==0)
            os << "const " << type;
          else 
            os << ", const " << type;
          if (!is_scalar) os << "*";
          os << " " << arg;
      }
      
      // Metadata arguments
      if (spec.exasim == true)
        os << ", const int modelnumber, const int N";
      else
        os << ", const int N, const int szf";            
      for (const std::string& arg : func.args) {
          if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
              os << ", const int sz" << arg;
          }
      }
      os << ")"<<"\"";
      if (i < spec.functions.size() - 1) os << ", ";
      os << "\n";
    }
    os << "    };\n\n";
    
    os << "    inputvectors = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";        

        // Collect vector arguments only (exclude scalars)
        std::vector<std::string> vector_args;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end()) {
                vector_args.push_back(arg);
            }
        }

        for (size_t i = 0; i < vector_args.size(); ++i) {
            const auto& name = vector_args[i];
            os << "{\"" << name << "\", " << name << "}";
            if (i + 1 < vector_args.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";

    os << "    inputscalars = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";        

        // Collect scalar arguments 
        std::vector<std::string> scalar_args;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) != spec.scalars.end()) {
                scalar_args.push_back(arg);
            }
        }

        for (size_t i = 0; i < scalar_args.size(); ++i) {
            const auto& name = scalar_args[i];
            os << "{\"" << name << "\", " << name << "}";
            if (i + 1 < scalar_args.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    jacobianInputs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";

        // Collect vector arguments that are in the Jacobian list
        std::vector<std::string> jacvecs;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end() &&
                std::find(spec.jacobian.begin(), spec.jacobian.end(), arg) != spec.jacobian.end()) {
                jacvecs.push_back(arg);
            }
        }

        for (size_t i = 0; i < jacvecs.size(); ++i) {
            os << jacvecs[i];
            if (i + 1 < jacvecs.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
    
    os << "    hessianInputs = {\n";
    for (size_t fi = 0; fi < spec.functions.size(); ++fi) {
        const auto& func = spec.functions[fi];
        os << "        {";

        // Collect vector arguments that are in the Jacobian list
        std::vector<std::string> jacvecs;
        for (const auto& arg : func.args) {
            if (std::find(spec.scalars.begin(), spec.scalars.end(), arg) == spec.scalars.end() &&
                std::find(spec.hessian.begin(), spec.hessian.end(), arg) != spec.hessian.end()) {
                jacvecs.push_back(arg);
            }
        }

        for (size_t i = 0; i < jacvecs.size(); ++i) {
            os << jacvecs[i];
            if (i + 1 < jacvecs.size()) os << ", ";
        }

        os << "}";
        if (fi + 1 < spec.functions.size()) os << ",";
        os << "\n";
    }
    os << "    };\n\n";
            
    os << "}\n\n";
}

void emitevaluateSymbolicFunctions(std::ostream& os, const ParsedSpec& spec) {
    os << "std::vector<Expression> SymbolicScalarsVectors::evaluateSymbolicFunctions(int call)\n";
    os << "{\n";
    os << "  std::vector<Expression> f;\n\n";
    os << "  switch (call) {\n";

    for (size_t i = 0; i < spec.functions.size(); ++i) {
        const auto& sig = spec.functions[i];
        os << "    case " << i << ":\n";
        os << "      f = " << sig.name << "(";
        for (size_t j = 0; j < sig.args.size(); ++j) {
            os << sig.args[j];
            if (j < sig.args.size() - 1)
                os << ", ";
        }
        os << ");\n";
        os << "      break;\n";
    }

    os << "    default:\n";
    os << "      throw std::runtime_error(\"Invalid function call in evaluateSymbolicFunctions\");\n";
    os << "  }\n\n";
    os << "  return f;\n";
    os << "}\n\n";
}
            
void emitfunc2cse(std::ostream& os) {
    os<<"void SymbolicScalarsVectors::func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f) {\n\n"
        "   vec_basic exprs;\n"
        "   for (const auto &fi : f) {\n"
        "       exprs.push_back(fi.get_basic());\n"
        "   }\n"
        "   cse(replacements, reduced_exprs, exprs);\n"
        "}\n\n";
}

void emitfuncjac2cse(std::ostream& os) {
    os << "void SymbolicScalarsVectors::funcjac2cse(vec_pair &replacements,\n"
       << "                                         vec_basic &reduced_exprs_f,\n"
       << "                                         std::vector<vec_basic> &reduced_exprs_J,\n"
       << "                                         const std::vector<Expression> &f,\n"
       << "                                         const std::vector<std::vector<Expression>>& inputs_J) {\n"
       << "    vec_basic exprs;\n\n"
       << "    // Track original sizes\n"
       << "    int n_f = static_cast<int>(f.size());\n"
       << "    std::vector<int> jacobian_sizes;\n\n"
       << "    // Add original function expressions\n"
       << "    for (const auto &fi : f) {\n"
       << "        exprs.push_back(fi.get_basic());\n"
       << "    }\n\n"
       << "    // Append Jacobians and record sizes for each input group\n"
       << "    for (const auto& input_vec : inputs_J) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        jacobian_sizes.push_back(m * n);\n"
       << "        for (int j = 0; j < n; ++j) {\n"
       << "            for (int i = 0; i < m; ++i) {\n"
       << "                exprs.push_back(f[i].diff(input_vec[j]).get_basic());\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"
       << "    // Apply CSE\n"
       << "    vec_basic reduced_exprs;\n"
       << "    cse(replacements, reduced_exprs, exprs);\n\n"
       << "    // Decompose: f\n"
       << "    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + n_f);\n\n"
       << "    // Decompose: Jacobian blocks\n"
       << "    int offset = n_f;\n"
       << "    for (const int sz : jacobian_sizes) {\n"
       << "        vec_basic block(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_J.push_back(block);\n"
       << "        offset += sz;\n"
       << "    }\n"
       << "}\n\n";
}

void emitfuncjachess2cse(std::ostream& os) {
    os << "void SymbolicScalarsVectors::funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,\n"
       << "     std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,\n"
       << "     const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,\n"
       << "     const std::vector<std::vector<Expression>>& inputs_H) {\n\n"

       << "    vec_basic exprs;\n"
       << "    int count_f = f.size();\n"
       << "    int count_J = 0;\n"
       << "    int count_H = 0;\n\n"

       << "    // Add original function expressions\n"
       << "    for (const auto &fi : f) {\n"
       << "        exprs.push_back(fi.get_basic());\n"
       << "    }\n\n"

       << "    // Compute and append all Jacobian entries for each input group\n"
       << "    std::vector<int> J_sizes;\n"
       << "    for (const auto& input_vec : inputs_J) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        J_sizes.push_back(m * n);\n"
       << "        count_J += m * n;\n"
       << "        for (int j = 0; j < n; ++j) {\n"
       << "            for (int i = 0; i < m; ++i) {\n"
       << "                exprs.push_back(f[i].diff(input_vec[j]).get_basic());\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"

       << "    // Compute and append all Hessian entries for each input group\n"
       << "    std::vector<int> H_sizes;\n"
       << "    for (const auto& input_vec : inputs_H) {\n"
       << "        int m = f.size();\n"
       << "        int n = input_vec.size();\n"
       << "        H_sizes.push_back(m * n * n);\n"
       << "        count_H += m * n * n;\n"
       << "        for (int k = 0; k < n; ++k) {\n"
       << "            for (int j = 0; j < n; ++j) {\n"
       //<< "                Expression first = f[i].diff(input_vec[j]);\n"
       << "                for (int i = 0; i < m; ++i) {\n"
       << "                    exprs.push_back(f[i].diff(input_vec[j]).diff(input_vec[k]).get_basic());\n"
       << "                }\n"
       << "            }\n"
       << "        }\n"
       << "    }\n\n"

       << "    // Apply Common Subexpression Elimination\n"
       << "    vec_basic reduced_exprs;\n"
       << "    cse(replacements, reduced_exprs, exprs);\n\n"

       << "    // Decompose reduced_exprs into f, J, and H\n"
       << "    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + count_f);\n\n"

       << "    int offset = count_f;\n"
       << "    for (int sz : J_sizes) {\n"
       << "        vec_basic Jblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_J.push_back(Jblock);\n"
       << "        offset += sz;\n"
       << "    }\n\n"

       << "    for (int sz : H_sizes) {\n"
       << "        vec_basic Hblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);\n"
       << "        reduced_exprs_H.push_back(Hblock);\n"
       << "        offset += sz;\n"
       << "    }\n"
       << "}\n\n";
}

void forloopstart(std::ostream& os, std::string framework)
{       
    if ((framework == "kokkos") || (framework == "Kokkos") || (framework == "KOKKOS")) {    
      os << "       cppfile << \"  Kokkos::parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, KOKKOS_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {    
      os << "       cppfile << \"  parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, GPU_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {    
      os << "       cppfile << \"  parallel_for(\\\"\"<< funcnames[functionid] <<\"\\\", N, GPU_LAMBDA(const size_t i) {\\n\";\n";                  
    }        
    else {
      os << "       cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n";   
    }
}

void forloopend(std::ostream& os, std::string framework)
{       
    if ((framework == "kokkos") || (framework == "Kokkos") || (framework == "KOKKOS")) {    
      os << "       cppfile << \"  });\\n\";\n";
    }
    else if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {       
      os << "       cppfile << \"  });\\n\";\n";
    }
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {        
      os << "       cppfile << \"  });\\n\";\n";
    }
    else {
      os << "       cppfile << \"  }\\n\";\n";   
    }
}

void emitCudaHipHpp(std::ostream& os, std::string framework)
{
    if ((framework == "cuda") || (framework == "Cuda") || (framework == "CUDA")) {             
      os << "    hfile << \"#include \\\"CudaHip.hpp\\\" \\n\\n\";\n";
    }
    else if ((framework == "hip") || (framework == "Hip") || (framework == "HIP")) {        
      os << "    hfile << \"#include \\\"CudaHip.hpp\\\" \\n\\n\";\n";
    }  
}

void emitfunc2cppfiles(std::ostream& os, const ParsedSpec& spec) {
  
    os << "void SymbolicScalarsVectors::func2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {\n";
        
    os << "    std::ios_base::openmode mode = std::ios::out;\n";
    os << "    if (append)\n";
    os << "        mode |= std::ios::app;\n";
    os << "    else\n";
    os << "        mode |= std::ios::trunc;\n\n";
    
    if (spec.exasim == false) {     
        os << "    // Generate function prototype header based on functionid\n";
        os << "    std::ofstream hfile(filename + std::string(\".h\"), mode);\n";
        os << "    hfile << \"#pragma once\\n\\n\";\n";
        emitCudaHipHpp(os, spec.framework);
        os << "    hfile << funcdecls[functionid] << \";\\n\";\n";
        os << "    hfile.close();\n";
        os << "\n\n";

        os << "    // Generate function prototype source based on functionid\n";
        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
        os << "    cppfile << \"#include \\\"\" << filename << \".h\\\"\\n\";\n";
        os << "    cppfile << \"#include <cmath>\\n\\n\";\n";
        os << "    cppfile << funcdecls[functionid] << \"\\n\";\n";    
    } else {
        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
        os << "    cppfile << \"void \" <<funcname << \"(" << spec.datatype << "* f, \";\n";
        os << "    cppfile << funcjacdecls[functionid] << \"\\n\";\n";      
    }
    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";           
    
    os << "   if (f.size() > 0) {\n";    

    os << "       vec_pair replacements;\n";
    os << "       vec_basic reduced_exprs;\n";
    os << "       func2cse(replacements, reduced_exprs, f);\n\n";

    os << "       // Determine variable usage\n";
    os << "       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "       for (const auto &expr : f) {\n";
    os << "           auto symbols = free_symbols(*expr.get_basic());\n";
    os << "           used.insert(symbols.begin(), symbols.end());\n";
    os << "       }\n\n";
    
    os << "       auto depends_on = [&](const Expression &sym) {\n";
    os << "           return used.count(sym.get_basic()) > 0;\n";
    os << "       };\n\n";    
    
    os << "       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "       C99CodePrinter cpp;\n";
            
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n";
    forloopstart(os, spec.framework);    
    os << "       \n";
    os << "       // Emit symbolic variable loads\n";
    os << "       for (const auto &[name, vec] : inputs) {\n";
    os << "           for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "               if (depends_on(vec[j])) { \n";
    os << "                 if (std::find(batch.begin(), batch.end(), name) != batch.end())\n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "                 else \n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"];\\n\";\n";
    os << "               }\n";
    os << "           }\n";
    os << "       }\n";
    os << "       \n";
    os << "       cppfile << \"\\n\";\n";
    os << "       \n";
    
//     for (size_t i = 0; i < sf.substitutions.size(); ++i) {
//         std::string var_name = cpp.apply(*sf.substitutions[i].first);
//         std::string rhs = cpp.apply(*sf.substitutions[i].second);
//         cppfile << "    double " << var_name << " = " << rhs << ";\n";
//     }
// 
                
    os << "       // Emit intermediate CSE substitutions\n";
    os << "       for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "           std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "           std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "           cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    //os << "        cppfile << \"    " << spec.datatype << " tmp\" << n << \" = \" << rhs << \";\\n\";\n";
    os << "       }\n";
    os << "       cppfile << \"\\n\";\n";
    os << "       \n";
    
    os << "       for (size_t n = 0; n < f.size(); ++n) {\n";
    //os << "        auto replaced = reduced_exprs[n]->subs(rename_map);\n";
    os << "           cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs[n]) << \";\\n\";\n";
    os << "       }\n";
    os << "       \n";
    //os << "    cppfile << \"  }\\n\";\n";  // forloop
    forloopend(os, spec.framework);
    os << "   }\n";  // os << "    if (f.size() > 0) {\n";         
    
    os << "    cppfile << \"}\\n\\n\";\n"; // end of cppfile function 
    os << "    cppfile.close();\n";
    os << "}\n\n";
    
//     for (int i = 0; i < k; ++i) {
//         cppfile << "    fb[" << i << "] = " << cpp.apply(sf.f[i].get_basic()) << ";\n";
//         for (int j = 0; j < m; ++j)
//             cppfile << "    dfdx_b[" << i*m + j << "] = " << cpp.apply(sf.dfdx[i][j].get_basic()) << ";\n";
//         for (int j = 0; j < n; ++j)
//             cppfile << "    dfdy_b[" << i*n + j << "] = " << cpp.apply(sf.dfdy[i][j].get_basic()) << ";\n";
//     }
    
//     os << "    // Emit final expressions\n";
//     os << "    map_basic_basic rename_map;\n";
//     os << "    for (size_t n = 0; n < replacements.size(); ++n)\n";
//     os << "        rename_map[replacements[n].first] = symbol(\"tmp\" + std::to_string(n));\n";
//     os << "    \n";
//     os << "    for (size_t n = 0; n < f.size(); ++n) {\n";
//     os << "        auto replaced = reduced_exprs[n]->subs(rename_map);\n";
//     os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "    }\n";
//     os << "    \n";
//     os << "    cppfile << \"  }\\n\";\n";
//     os << "    cppfile << \"}\\n\\n\";\n";    
//     os << "}\n\n";
}

void emitfuncjac2cppfiles(std::ostream& os, const ParsedSpec& spec) {
    os << "void SymbolicScalarsVectors::funcjac2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {\n";
    
    os << "    std::ios_base::openmode mode = std::ios::out;\n";
    os << "    if (append)\n";
    os << "        mode |= std::ios::app;\n";
    os << "    else\n";
    os << "        mode |= std::ios::trunc;\n\n";
    
    if (spec.exasim == false) {     
        os << "    // Generate function prototype header based on functionid\n";
        os << "    std::ofstream hfile(filename + std::string(\".h\"), mode);\n";
        os << "    hfile << \"void \" << funcname << \"jac\" << \"(" << spec.datatype << "* f, \";\n";
        os << "    int nJ = jacobianInputs[functionid].size();\n";
        os << "    for (int k = 0; k < nJ; ++k)\n";
        os << "        hfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
        os << "    hfile << funcjacdecls[functionid] << \";\\n\";\n";    
        os << "    hfile.close();\n";
        os << "\n";    

        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
        os << "    cppfile << \"void \" << funcname << \"jac\" << \"(" << spec.datatype << "* f, \";\n";
    } else {
        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";   
        os << "    cppfile << \"void \" << funcname << \"(" << spec.datatype << "* f, \";\n";        
        os << "    int nJ = jacobianInputs[functionid].size();\n";
        os << "    for (int k = 0; k < nJ; ++k)\n";
        os << "        cppfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
        os << "    cppfile << funcjacdecls[functionid] << \"\\n\";\n";        
    }
    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";        
    
    os << "   if (f.size() > 0) {\n";      
    os << "       vec_pair replacements;\n";
    os << "       vec_basic reduced_exprs_f;\n";
    os << "       std::vector<vec_basic> reduced_exprs_J;\n";
    os << "       funcjac2cse(replacements, reduced_exprs_f, reduced_exprs_J, f, jacobianInputs[functionid]);\n\n";
    
    os << "       // Determine variable usage\n";
    os << "       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "       for (const auto &expr : f) {\n";
    os << "           auto symbols = free_symbols(*expr.get_basic());\n";
    os << "           used.insert(symbols.begin(), symbols.end());\n";
    os << "       }\n\n";
    
    os << "       auto depends_on = [&](const Expression &sym) {\n";
    os << "           return used.count(sym.get_basic()) > 0;\n";
    os << "       };\n\n";
    
    os << "       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "       C99CodePrinter cpp;\n";    
    
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n\n";
    forloopstart(os, spec.framework);    
    os << "       for (const auto &[name, vec] : inputs) {\n";
    os << "           for (size_t j = 0; j < vec.size(); ++j) {\n";    
    os << "               if (depends_on(vec[j])) { \n";
    os << "                 if (std::find(batch.begin(), batch.end(), name) != batch.end())\n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "                 else \n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"];\\n\";\n";
    os << "               }\n";        
    os << "           }\n";
    os << "       }\n\n";

    os << "       for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "           std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "           std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "           cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    os << "       }\n";
    os << "       cppfile << \"\\n\";\n";
    os << "       \n";

    // No rename_map or substitution
    os << "       for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
    os << "           cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs_f[n]) << \";\\n\";\n";
    os << "       }\n\n";

    os << "       for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
    os << "           for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
    os << "               cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_J[k][j]) << \";\\n\";\n";
    os << "           }\n";
    os << "       }\n";    
    //os << "    cppfile << \"  }\\n\";\n";
    forloopend(os, spec.framework);
    os << "   }\n";  // os << "    if (f.size() > 0) {\n";        
    
    os << "    cppfile << \"}\\n\\n\";\n";
    os << "    cppfile.close();\n";
    os << "}\n\n";
}

void emitfuncjachess2cppfiles(std::ostream& os, const ParsedSpec& spec) {
    os << "void SymbolicScalarsVectors::funcjachess2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {\n";
    os << "    vec_pair replacements;\n";
    os << "    vec_basic reduced_exprs_f;\n";
    os << "    std::vector<vec_basic> reduced_exprs_J, reduced_exprs_H;\n";
    os << "    funcjachess2cse(replacements, reduced_exprs_f, reduced_exprs_J, reduced_exprs_H, f, jacobianInputs[functionid], hessianInputs[functionid]);\n\n";

    os << "    // Determine variable usage\n";
    os << "    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "    for (const auto &expr : f) {\n";
    os << "        auto symbols = free_symbols(*expr.get_basic());\n";
    os << "        used.insert(symbols.begin(), symbols.end());\n";
    os << "    }\n\n";
    
    os << "    auto depends_on = [&](const Expression &sym) {\n";
    os << "        return used.count(sym.get_basic()) > 0;\n";
    os << "    };\n\n";

    os << "    std::ios_base::openmode mode = std::ios::out;\n";
    os << "    if (append)\n";
    os << "        mode |= std::ios::app;\n";
    os << "    else\n";
    os << "        mode |= std::ios::trunc;\n\n";
    
    os << "    // Generate function prototype header based on functionid\n";
    os << "    std::ofstream hfile(filename + std::string(\".h\"), mode);\n";    
    os << "    hfile << \"void \" << funcname << \"jachess\" << \"(" << spec.datatype << "* f, \";\n";
    os << "    int nJ = reduced_exprs_J.size();\n";
    os << "    int nH = reduced_exprs_H.size();\n";
    os << "    for (int k = 0; k < nJ; ++k) hfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    for (int k = 0; k < nH; ++k) hfile << \""<< spec.datatype <<"* H\" << (k+1) << \", \";\n";    
    os << "    hfile << funcjacdecls[functionid] << \";\\n\";\n";    
    os << "    hfile.close();\n";
    os << "\n";    
            
    os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
    os << "    cppfile << \"void \" << funcname << \"jachess\" << \"(" << spec.datatype << "* f, \";\n";    
    os << "    for (int k = 0; k < nJ; ++k) cppfile << \""<< spec.datatype <<"* J\" << (k+1) << \", \";\n";    
    os << "    for (int k = 0; k < nH; ++k) cppfile << \""<< spec.datatype <<"* H\" << (k+1) << \", \";\n";    
    os << "    cppfile << funcjacdecls[functionid] << \"\\n\";\n";    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";

//     os << "    const std::vector<std::string>& args = funcargs[functionid];\n";
//     os << "    const std::vector<std::string>& sizes = funcargssizes[functionid];\n";
//     os << "    for (size_t i = 0; i < args.size(); ++i)\n";
//     os << "        cppfile << \"const double* \" << args[i] << \", \";\n";
//     os << "    cppfile << \"const int N, const int szf, \";\n";
//     os << "    for (size_t i = 0; i < sizes.size(); ++i) {\n";
//     os << "        cppfile << \"const int \" << sizes[i];\n";
//     os << "        if (i + 1 < sizes.size()) cppfile << \", \";\n";
//     os << "    }\n";
//     os << "    cppfile << \")\\n\";\n";
//     os << "    cppfile << \"{\\n\\n\";\n\n";

    os << "    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "    C99CodePrinter cpp;\n";    
    
    //os << "    cppfile << \"  for (int i = 0; i < N; ++i) {\\n\";\n\n";
    forloopstart(os, spec.framework);    
    os << "    for (const auto &[name, vec] : inputs) {\n";
    os << "        for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "            if (depends_on(vec[j])) { \n";
    os << "              if (std::find(batch.begin(), batch.end(), name) != batch.end())\n";
    os << "                 cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "               else \n";
    os << "                 cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"];\\n\";\n";
    os << "            }\n";
    os << "        }\n";
    os << "    }\n\n";

//     os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
//     os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
//     os << "        cppfile << \"    " << spec.datatype << " tmp\" << n << \" = \" << rhs << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    map_basic_basic rename_map;\n";
//     os << "    for (size_t n = 0; n < replacements.size(); ++n)\n";
//     os << "        rename_map[replacements[n].first] = symbol(\"tmp\" + std::to_string(n));\n\n";
// 
//     os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
//     os << "        auto replaced = reduced_exprs_f[n]->subs(rename_map);\n";
//     os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "    }\n\n";
// 
//     os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
//     os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
//     os << "            auto replaced = reduced_exprs_J[k][j]->subs(rename_map);\n";
//     os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "        }\n";
//     os << "    }\n\n";
// 
//     os << "    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {\n";
//     os << "        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {\n";
//     os << "            auto replaced = reduced_exprs_H[k][j]->subs(rename_map);\n";
//     os << "            cppfile << \"    H\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*replaced) << \";\\n\";\n";
//     os << "        }\n";
//     os << "    }\n";

    os << "    for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "        std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "        std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "        cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    os << "    }\n";
    os << "    cppfile << \"\\n\";\n";
    os << "    \n";

    // Skip substitution, directly print the reduced expressions
    os << "    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {\n";
    os << "        cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs_f[n]) << \";\\n\";\n";
    os << "    }\n\n";

    os << "    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {\n";
    os << "        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {\n";
    os << "            cppfile << \"    J\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_J[k][j]) << \";\\n\";\n";
    os << "        }\n";
    os << "    }\n\n";

    os << "    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {\n";
    os << "        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {\n";
    os << "            cppfile << \"    H\" << (k+1) << \"[\" << j << \" * N + i] = \" << cpp.apply(*reduced_exprs_H[k][j]) << \";\\n\";\n";
    os << "        }\n";
    os << "    }\n";    
    //os << "    cppfile << \"  }\\n\";\n";
    forloopend(os, spec.framework);
    
    os << "    cppfile << \"}\\n\\n\";\n";
    os << "    cppfile.close();\n";
    os << "}\n\n";
}


void emitinitfunc2cppfiles(std::ostream& os, const ParsedSpec& spec) {
  
    os << "void SymbolicScalarsVectors::initfunc2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append, int framework) {\n";
        
    os << "    std::ios_base::openmode mode = std::ios::out;\n";
    os << "    if (append)\n";
    os << "        mode |= std::ios::app;\n";
    os << "    else\n";
    os << "        mode |= std::ios::trunc;\n\n";
    
    if (spec.exasim == false) {     
        os << "    // Generate function prototype header based on functionid\n";
        os << "    std::ofstream hfile(filename + std::string(\".h\"), mode);\n";
        os << "    hfile << \"#pragma once\\n\\n\";\n";
        emitCudaHipHpp(os, spec.framework);
        os << "    hfile << funcdecls[functionid] << \";\\n\";\n";
        os << "    hfile.close();\n";
        os << "\n\n";

        os << "    // Generate function prototype source based on functionid\n";
        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
        os << "    cppfile << \"#include \\\"\" << filename << \".h\\\"\\n\";\n";
        os << "    cppfile << \"#include <cmath>\\n\\n\";\n";
        os << "    cppfile << funcdecls[functionid] << \"\\n\";\n";    
    } else {
        os << "    std::ofstream cppfile(filename + std::string(\".cpp\"), mode);\n";
        os << "    cppfile << \"void \" <<funcname << \"(" << spec.datatype << "* f, \";\n";
        os << "    cppfile << \"const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)\" << \"\\n\";\n";      
    }
    
    os << "    cppfile << \"{\\n\\n\";";
    os << "\n\n";           
    
    os << "   if (f.size() > 0) {\n";    

    os << "       vec_pair replacements;\n";
    os << "       vec_basic reduced_exprs;\n";
    os << "       func2cse(replacements, reduced_exprs, f);\n\n";

    os << "       // Determine variable usage\n";
    os << "       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;\n";
    os << "       for (const auto &expr : f) {\n";
    os << "           auto symbols = free_symbols(*expr.get_basic());\n";
    os << "           used.insert(symbols.begin(), symbols.end());\n";
    os << "       }\n\n";
    
    os << "       auto depends_on = [&](const Expression &sym) {\n";
    os << "           return used.count(sym.get_basic()) > 0;\n";
    os << "       };\n\n";    
    
    os << "       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];\n";
    os << "       C99CodePrinter cpp;\n";
            
    os << "       if (framework==0) \n";    
    forloopstart(os, "C++");    
    os << "       else if (framework==1) \n";     
    forloopstart(os, "Kokkos");    
    os << "       else if (framework==2) \n";      
    forloopstart(os, "Cuda");    
    os << "       cppfile << \"    int p = i%npe; \\n\";\n";   
    os << "       cppfile << \"    int e = i/npe; \\n\";\n";   
    os << "       \n";
    os << "       // Emit symbolic variable loads\n";
    os << "       for (const auto &[name, vec] : inputs) {\n";
    os << "           for (size_t j = 0; j < vec.size(); ++j) {\n";
    os << "               if (depends_on(vec[j])) { \n";
    os << "                 if (std::find(batch.begin(), batch.end(), name) != batch.end())\n";
    //os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*N+i];\\n\";\n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"*npe+p+npe*ncx*e];\\n\";\n";
    os << "                 else \n";
    os << "                   cppfile << \"    " << spec.datatype << " \" << name << j << \" = \" << name << \"[\" << j << \"];\\n\";\n";
    os << "               }\n";    
    os << "           }\n";
    os << "       }\n";
    os << "       \n";
    os << "       cppfile << \"\\n\";\n";
    os << "       \n";
                   
    os << "       // Emit intermediate CSE substitutions\n";
    os << "       for (size_t n = 0; n < replacements.size(); ++n) {\n";
    os << "           std::string var_name = cpp.apply(*replacements[n].first);\n";
    os << "           std::string rhs = cpp.apply(*replacements[n].second);\n";
    os << "           cppfile << \"    " << spec.datatype << " \" << var_name << \" = \" << rhs << \";\\n\";\n";
    os << "       }\n";
    os << "       cppfile << \"\\n\";\n";
    os << "       \n";
    
    os << "       for (size_t n = 0; n < f.size(); ++n) {\n";    
    os << "           cppfile << \"    f[p+npe*\" << n << \" +npe*nce*e] = \" << cpp.apply(*reduced_exprs[n]) << \";\\n\";\n";
    //os << "           cppfile << \"    f[\" << n << \" * N + i] = \" << cpp.apply(*reduced_exprs[n]) << \";\\n\";\n";
    os << "       }\n";
    os << "       \n";   
    os << "       if (framework==0) \n";      
    forloopend(os, "C++");
    os << "       else if (framework==1) \n";     
    forloopend(os, "Kokkos");
    os << "       else if (framework==2) \n";      
    forloopend(os, "Cuda");
    os << "   }\n";  // os << "    if (f.size() > 0) {\n";         
    
    os << "    cppfile << \"}\\n\\n\";\n"; // end of cppfile function 
    os << "    cppfile.close();\n";
    os << "}\n\n";    
}

void emitUbouFbou(std::ostream& os) {
    os << "void SymbolicScalarsVectors::appendUbouFbou(const std::string& filename, const std::string& funcname, int nbc) {\n";
    os << "    std::ostringstream tmp;\n\n";
    os << "    tmp << \"void \" << funcname << \"(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\\n\";\n";
    os << "        tmp << \"           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\\n\";\n";
    os << "        tmp << \"           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\\n\";\n";
    os << "        tmp << \"           const int ncx, const int nco, const int ncw) {\\n\";\n\n";

    os << "    for (int k = 1; k <= nbc; ++k) {\n";
    os << "        if (k == 1)\n";
    os << "            tmp << \"    if (ib == 1 )\\n\";\n";
    os << "        else\n";
    os << "            tmp << \"    else if (ib == \" << k << \" )\\n\";\n";
    os << "        tmp << \"        \" << funcname << k << \"(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,\\n\";\n";
    os << "        tmp << \"                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);\\n\";\n";
    os << "    }\n\n";
    
    os << "    tmp << \"}\\n\";\n\n";
    
    os << "    std::ofstream cppfile(filename + \".cpp\", std::ios::out | std::ios::app);\n";    
    os << "    cppfile << tmp.str();\n";
    os << "    cppfile.close();\n";
    os << "}\n";
}

void emitFbouHdg(std::ostream& os) {
    os << "void SymbolicScalarsVectors::appendFbouHdg(const std::string& filename, const std::string& funcname, int nbc) {\n";
    os << "    std::ostringstream tmp;\n\n";
    os << "    tmp << \"void \" << funcname << \"(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\\n\";\n";
    os << "        tmp << \"           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\\n\";\n";
    os << "        tmp << \"           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\\n\";\n";
    os << "        tmp << \"           const int ncx, const int nco, const int ncw) {\\n\";\n\n";

    os << "    for (int k = 1; k <= nbc; ++k) {\n";
    os << "        if (k == 1)\n";
    os << "            tmp << \"    if (ib == 1 )\\n\";\n";
    os << "        else\n";
    os << "            tmp << \"    else if (ib == \" << k << \" )\\n\";\n";
    os << "        tmp << \"        \" << funcname << k << \"(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,\\n\";\n";
    os << "        tmp << \"                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);\\n\";\n";
    os << "    }\n\n";
    
    os << "    tmp << \"}\\n\";\n\n";
    
    os << "    std::ofstream cppfile(filename + \".cpp\", std::ios::out | std::ios::app);\n";    
    os << "    cppfile << tmp.str();\n";
    os << "    cppfile.close();\n";
    os << "}\n";
}

void CodeGenerator::generateSymbolicScalarsVectorsCpp(const std::string& filename) const {
    std::ofstream os(filename);

    os << "#include \"SymbolicScalarsVectors.hpp\"\n\n";
        
    emitSymbolicScalarsVectors(os, spec);
    
    emitevaluateSymbolicFunctions(os, spec);
                           
    emitfunc2cse(os);
        
    emitfuncjac2cse(os);
    
    emitfuncjachess2cse(os);    
       
    emitfunc2cppfiles(os, spec);    
    
    emitfuncjac2cppfiles(os, spec);       
    
    emitfuncjachess2cppfiles(os, spec);
    
    emitinitfunc2cppfiles(os, spec);    
    
    emitUbouFbou(os);
    
    emitFbouHdg(os);
    
    os.close();  
}

void CodeGenerator::generateEmptySourcewCpp(std::string modelpath) const {  
  {
    std::ofstream os(make_path(modelpath,  "KokkosSourcew.cpp"));
    os << "void KokkosSourcew(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();
  }
  {
    std::ofstream os(make_path(modelpath,  "HdgSourcewonly.cpp"));
    os << "void HdgSourcewonly(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();    
  }
  
    std::ofstream os2(make_path(modelpath,  "HdgSourcew.cpp"));
    os2 << "void HdgSourcew(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os2 << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os2 << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();       
}

void CodeGenerator::generateEmptyOutputCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosOutput.cpp"));
    os << "void KokkosOutput(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyMonitorCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosMonitor.cpp"));
    os << "void KokkosMonitor(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyAvfieldCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosAvfield.cpp"));
    os << "void KokkosAvfield(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyInituCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosInitu.cpp"));
    os << "void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
    
    std::ofstream os2(make_path(modelpath,  "cpuInitu.cpp"));
    os2 << "void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyInitqCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosInitq.cpp"));
    os << "void KokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
    
    std::ofstream os2(make_path(modelpath,  "cpuInitq.cpp"));
    os2 << "void cpuInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyInitudgCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosInitudg.cpp"));
    os << "void KokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
    
    std::ofstream os2(make_path(modelpath,  "cpuInitudg.cpp"));
    os2 << "void cpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyInitodgCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosInitodg.cpp"));
    os << "void KokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
    
    std::ofstream os2(make_path(modelpath,  "cpuInitodg.cpp"));
    os2 << "void cpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyInitwdgCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosInitwdg.cpp"));
    os << "void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
    
    std::ofstream os2(make_path(modelpath,  "cpuInitwdg.cpp"));
    os2 << "void cpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyFintCpp(std::string modelpath) const {  
  {
    std::ofstream os(make_path(modelpath,  "KokkosFint.cpp"));
    os << "void KokkosFint(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "             const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
    os << "             const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx,\n";
    os << "             const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
  }
    
  {
    std::ofstream os(make_path(modelpath,  "HdgFintonly.cpp"));
    os << "void HdgFintonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "             const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
    os << "             const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx,\n";
    os << "             const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
  }
  
    std::ofstream os2(make_path(modelpath,  "HdgFint.cpp"));
    os2 << "void HdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os2 << "             const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
    os2 << "             const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx,\n";
    os2 << "             const int nco, const int ncw)\n";
    os2 << "{\n";
    os2 << "}\n";
    os2.close();        
}

void CodeGenerator::generateEmptyFhatCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosFhat.cpp"));
    os << "void KokkosFhat(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyUhatCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosUhat.cpp"));
    os << "void KokkosUhat(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyStabCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosStab.cpp"));
    os << "void KokkosStab(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();        
}

void CodeGenerator::generateEmptyEoSCpp(std::string modelpath) const {  
    // KokkosEoS.cpp
    {
        std::ofstream os(make_path(modelpath,  "KokkosEoS.cpp"));
        os << "void KokkosEoS(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, "
           << "const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, "
           << "const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, "
           << "const int nce, const int npe, const int ne)\n"
           << "{\n"
           << "}\n";
    }

    // KokkosEoSdu.cpp
    {
        std::ofstream os(make_path(modelpath,  "KokkosEoSdu.cpp"));
        os << "void KokkosEoSdu(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, "
           << "const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, "
           << "const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, "
           << "const int nce, const int npe, const int ne)\n"
           << "{\n"
           << "}\n";
    }

    // KokkosEoSdw.cpp
    {
        std::ofstream os(make_path(modelpath,  "KokkosEoSdw.cpp"));
        os << "void KokkosEoSdw(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, "
           << "const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, "
           << "const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, "
           << "const int nce, const int npe, const int ne)\n"
           << "{\n"
           << "}\n";
    }

    // HdgEoS.cpp
    {
        std::ofstream os(make_path(modelpath,  "HdgEoS.cpp"));
        os << "void HdgEoS(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, "
           << "const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, "
           << "const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, "
           << "const int nco, const int ncw)\n"
           << "{\n"
           << "}\n";
    }
}

void CodeGenerator::generateEmptyVisScalarsCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosVisScalars.cpp"));
    os << "void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();      
}

void CodeGenerator::generateEmptyVisVectorsCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosVisVectors.cpp"));
    os << "void KokkosVisVectors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();      
}

void CodeGenerator::generateEmptyVisTensorsCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosVisTensors.cpp"));
    os << "void KokkosVisTensors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();      
}

void CodeGenerator::generateEmptyQoIvolumeCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "KokkosQoIvolume.cpp"));
    os << "void KokkosQoIvolume(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "                 const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc,\n";
    os << "                 const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();      
}

void CodeGenerator::generateEmptyQoIboundaryCpp(std::string modelpath) const {    
    std::ofstream os(make_path(modelpath,  "KokkosQoIboundary.cpp"));
    os << "void KokkosQoIboundary(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg,\n";
    os << "             const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
    os << "             const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx,\n";
    os << "             const int nco, const int ncw)\n";
    os << "{\n";
    os << "}\n";
    os.close();          
}

void CodeGenerator::generateLibPDEModelHpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath, "libpdemodel.hpp"));

    os << "#pragma once\n\n";

    os << "void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosAvfield(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosEoS(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosEoSdu(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosEoSdw(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosFhat(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void KokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void KokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void KokkosMonitor(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosOutput(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosSourcew(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne);\n";
    os << "void KokkosStab(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosTdfunc(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosUbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosUhat(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";

    // CPU initialization routines
    os << "void cpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void cpuInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void cpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";
    os << "void cpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne);\n";

    // HDG functions
    os << "void HdgEoS(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgFintonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgSourcew(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void HdgSourcewonly(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";

    os << "void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosVisVectors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosVisTensors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosQoIvolume(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";
    os << "void KokkosQoIboundary(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw);\n";

    os.close(); 
}

void CodeGenerator::generateLibPDEModelCpp(std::string modelpath) const {  
    std::ofstream os(make_path(modelpath,  "libpdemodel.cpp"));

    os << "#include <cmath>\n";
    os << "#include <Kokkos_Core.hpp>\n\n";
    
    os << "#ifdef USE_FLOAT\n";
    os << "typedef float dstype;\n";
    os << "#else\n";
    os << "typedef double dstype; //  double is default precision \n";
    os << "#endif\n\n";
        
    os << "using namespace std;\n\n";
    
    os << "#include \"KokkosFlux.cpp\"\n";
    os << "#include \"KokkosFhat.cpp\"\n";
    os << "#include \"KokkosFbou.cpp\"\n";
    os << "#include \"KokkosUbou.cpp\"\n";
    os << "#include \"KokkosUhat.cpp\"\n";
    os << "#include \"KokkosStab.cpp\"\n";
    os << "#include \"KokkosSource.cpp\"\n";
    os << "#include \"KokkosSourcew.cpp\"\n";
    os << "#include \"KokkosOutput.cpp\"\n";
    os << "#include \"KokkosMonitor.cpp\"\n";
    os << "#include \"KokkosInitu.cpp\"\n";
    os << "#include \"KokkosInitq.cpp\"\n";
    os << "#include \"KokkosInitwdg.cpp\"\n";
    os << "#include \"KokkosInitudg.cpp\"\n";
    os << "#include \"KokkosInitodg.cpp\"\n";
    os << "#include \"KokkosEoS.cpp\"\n";
    os << "#include \"KokkosEoSdu.cpp\"\n";
    os << "#include \"KokkosEoSdw.cpp\"\n";
    os << "#include \"KokkosAvfield.cpp\"\n";
    os << "#include \"KokkosTdfunc.cpp\"\n\n";

    os << "#include \"cpuInitu.cpp\"\n";
    os << "#include \"cpuInitq.cpp\"\n";
    os << "#include \"cpuInitwdg.cpp\"\n";
    os << "#include \"cpuInitudg.cpp\"\n";
    os << "#include \"cpuInitodg.cpp\"\n\n";

    os << "#include \"HdgFlux.cpp\"\n";
    os << "#include \"HdgFbou.cpp\"\n";
    os << "#include \"HdgFbouonly.cpp\"\n";
    os << "#include \"HdgFint.cpp\"\n";
    os << "#include \"HdgFintonly.cpp\"\n";
    os << "#include \"HdgSource.cpp\"\n";
    os << "#include \"HdgSourcew.cpp\"\n";
    os << "#include \"HdgSourcewonly.cpp\"\n";
    os << "#include \"HdgEoS.cpp\"\n";

    os << "#include \"KokkosVisScalars.cpp\"\n";
    os << "#include \"KokkosVisVectors.cpp\"\n";
    os << "#include \"KokkosVisTensors.cpp\"\n";
    os << "#include \"KokkosQoIvolume.cpp\"\n";
    os << "#include \"KokkosQoIboundary.cpp\"\n";

    os.close(); 
}
