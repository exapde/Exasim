#pragma once

#include "SymbolicFunctions.hpp"

class SymbolicScalarsVectors {

public:

    // input symbolic scalars
    Expression t;

    // input symbolic vectors
    std::vector<Expression> tau;
    std::vector<Expression> w;
    std::vector<Expression> eta;
    std::vector<Expression> uq;
    std::vector<Expression> n;
    std::vector<Expression> uhat;
    std::vector<Expression> v;
    std::vector<Expression> mu;
    std::vector<Expression> x;

    // vector sizes
    int sztau;
    int szw;
    int szeta;
    int szuq;
    int szn;
    int szuhat;
    int szv;
    int szmu;
    int szx;
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
};
