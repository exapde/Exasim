#pragma once

#include "SymbolicFunctions.hpp"
#include <algorithm>
#include <tuple>

class SymbolicScalarsVectors {

public:

    // path to model folder; Code2Cpp's main() sets this from argv[1].
    // Default is empty (writes relative to cwd) so the committed
    // snapshot doesn't bake in any one developer's absolute path.
    std::string modelpath = "";

    // input symbolic scalars
    Expression t;

    // input symbolic vectors
    std::vector<Expression> x;
    std::vector<Expression> mu;
    std::vector<Expression> v;
    std::vector<Expression> uhat;
    std::vector<Expression> n;
    std::vector<Expression> uq;
    std::vector<Expression> eta;
    std::vector<Expression> w;
    std::vector<Expression> tau;
    std::vector<Expression> uext;

    // vector sizes
    int szx;
    int szmu;
    int szv;
    int szuhat;
    int szn;
    int szuq;
    int szeta;
    int szw;
    int sztau;
    int szuext;
    bool exasim;

    std::vector<bool> outputfunctions;
    std::vector<std::vector<std::string>> funcargs;
    std::vector<std::vector<std::string>> funcargssizes;
    std::vector<std::string> funcnames;
    std::vector<std::string> funcdecls;
    std::vector<std::string> funcjacdecls;

    std::vector<std::vector<std::pair<std::string, std::vector<Expression>>>> inputvectors;
    std::vector<std::vector<std::pair<std::string, Expression>>> inputscalars;
    std::vector<std::vector<std::vector<Expression>>> jacobianInputs;
    std::vector<std::vector<std::vector<Expression>>> hessianInputs;

    std::vector<std::string> batch;
    SymbolicScalarsVectors();

    std::vector<Expression> evaluateSymbolicFunctions(int call);

    void func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f);

    void funcjac2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,
                         std::vector<vec_basic> &reduced_exprs_J, const std::vector<Expression> &f,
                         const std::vector<std::vector<Expression>>& inputs_J);

    void funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,
                         std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,
                         const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,
                         const std::vector<std::vector<Expression>>& inputs_H);

    void func2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);
    void funcjac2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);
    void funcjachess2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append);
    void initfunc2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append, int framework);
    void appendUbouFbou(const std::string& filename, const std::string& funcname, int nbc);
    void appendFbouHdg(const std::string& filename, const std::string& funcname, int nbc);
    void appendFextonly(const std::string& filename, const std::string& funcname, int nbc);
    void appendFext(const std::string& filename, const std::string& funcname, int nbc);

    // Emit a single header-only `my_model.hpp` consumable by
    // `<exasim/model.hpp>`'s templated FEM internals (HOT.4).
    void emit_pointwise_value(std::ostream& os, const std::string& method_name,
                              const std::string& cpp_signature,
                              const std::vector<Expression>& f, int functionid);
    void emit_pointwise_value_per_ib(std::ostream& os, const std::string& method_name,
                                     const std::string& cpp_signature,
                                     const std::vector<Expression>& f, int functionid, int szuhat);
    void emit_pointwise_value_and_jac(std::ostream& os, const std::string& value_method_name,
                                     const std::vector<std::string>& jac_method_names,
                                     const std::string& cpp_signature_value,
                                     const std::vector<std::string>& cpp_signatures_jac,
                                     const std::vector<Expression>& f, int functionid);
    void generateModelHeader(const std::string& filename);
};
