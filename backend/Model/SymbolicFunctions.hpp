#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/printers/codegen.h>
#include <symengine/symbol.h>
#include <symengine/matrix.h>

using SymEngine::Expression;
using SymEngine::vec_pair;
using SymEngine::vec_basic;
using SymEngine::symbol;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::map_basic_basic;
using SymEngine::CodePrinter;
using SymEngine::C99CodePrinter;

std::vector<Expression> Flux(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> Source(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> Tdfunc(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> Fbou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> Ubou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> FbouHdg(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> Initu(const std::vector<Expression>& x, const std::vector<Expression>& eta, const std::vector<Expression>& mu);
std::vector<Expression> VisScalars(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
std::vector<Expression> VisVectors(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t);
