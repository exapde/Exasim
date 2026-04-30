#include "SymbolicScalarsVectors.hpp"

SymbolicScalarsVectors::SymbolicScalarsVectors() {

    t =   Expression("t");

    szx = 2;
    x.resize(2);
    for (int i = 0; i < 2; ++i) {
         x[i] = Expression("x"  + std::to_string(i));
    }

    szmu = 1;
    mu.resize(1);
    for (int i = 0; i < 1; ++i) {
         mu[i] = Expression("mu"  + std::to_string(i));
    }

    szv = 0;
    v.resize(0);
    for (int i = 0; i < 0; ++i) {
         v[i] = Expression("v"  + std::to_string(i));
    }

    szuhat = 1;
    uhat.resize(1);
    for (int i = 0; i < 1; ++i) {
         uhat[i] = Expression("uhat"  + std::to_string(i));
    }

    szn = 2;
    n.resize(2);
    for (int i = 0; i < 2; ++i) {
         n[i] = Expression("n"  + std::to_string(i));
    }

    szuq = 3;
    uq.resize(3);
    for (int i = 0; i < 3; ++i) {
         uq[i] = Expression("uq"  + std::to_string(i));
    }

    szeta = 0;
    eta.resize(0);
    for (int i = 0; i < 0; ++i) {
         eta[i] = Expression("eta"  + std::to_string(i));
    }

    szw = 0;
    w.resize(0);
    for (int i = 0; i < 0; ++i) {
         w[i] = Expression("w"  + std::to_string(i));
    }

    sztau = 1;
    tau.resize(1);
    for (int i = 0; i < 1; ++i) {
         tau[i] = Expression("tau"  + std::to_string(i));
    }

    szuext = 1;
    uext.resize(1);
    for (int i = 0; i < 1; ++i) {
         uext[i] = Expression("uext"  + std::to_string(i));
    }

    exasim = true;

    outputfunctions.assign(13, false);
    outputfunctions[0] = true;
    outputfunctions[1] = true;
    outputfunctions[2] = true;
    outputfunctions[3] = true;
    outputfunctions[4] = true;
    outputfunctions[5] = true;
    outputfunctions[6] = true;
    outputfunctions[7] = true;
    outputfunctions[8] = true;
    outputfunctions[9] = true;
    outputfunctions[10] = true;
    outputfunctions[11] = true;
    outputfunctions[12] = true;

    batch = {"x", "uq", "v", "w", "uhat", "n", "uext"};

    funcnames = {"Flux", "Source", "Tdfunc", "Fbou", "Ubou", "FbouHdg", "Initu", "VisScalars", "VisVectors", "QoIvolume", "QoIboundary", "Fint", "Fext"};

    funcargs = {
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "tau", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "tau", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "tau", "eta", "mu", "t"},
        {"x", "eta", "mu"},
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "tau", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "tau", "eta", "mu", "t"},
        {"x", "uq", "v", "w", "uhat", "n", "uext", "tau", "eta", "mu", "t"}
    };

    funcargssizes = {
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "sztau", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "sztau", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "sztau", "szeta", "szmu", "szt"},
        {"szx", "szeta", "szmu"},
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "sztau", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "sztau", "szeta", "szmu", "szt"},
        {"szx", "szuq", "szv", "szw", "szuhat", "szn", "szuext", "sztau", "szeta", "szmu", "szt"}
    };

    funcdecls = {
       "void Flux(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void Source(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void Tdfunc(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void Fbou(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "void Ubou(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "void FbouHdg(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "void Initu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int szx, const int szeta, const int szmu)", 
       "void VisScalars(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void VisVectors(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void QoIvolume(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "void QoIboundary(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "void Fint(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "void Fext(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* uext, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int szuext, const int sztau, const int szeta, const int szmu)"
    };

    funcjacdecls = {
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int szx, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)", 
       "const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* uext, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int szuext, const int sztau, const int szeta, const int szmu)"
    };

    inputvectors = {
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"tau", tau}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"tau", tau}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"tau", tau}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"tau", tau}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"tau", tau}, {"eta", eta}, {"mu", mu}},
        {{"x", x}, {"uq", uq}, {"v", v}, {"w", w}, {"uhat", uhat}, {"n", n}, {"uext", uext}, {"tau", tau}, {"eta", eta}, {"mu", mu}}
    };

    inputscalars = {
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}},
        {{"t", t}}
    };

    jacobianInputs = {
        {uq, w},
        {uq, w},
        {uq, w},
        {uq, w, uhat},
        {uq, w, uhat},
        {uq, w, uhat},
        {},
        {uq, w},
        {uq, w},
        {uq, w},
        {uq, w, uhat},
        {uq, w, uhat},
        {uq, w, uhat}
    };

    hessianInputs = {
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {}
    };

}

std::vector<Expression> SymbolicScalarsVectors::evaluateSymbolicFunctions(int call)
{
  std::vector<Expression> f;

  switch (call) {
    case 0:
      f = Flux(x, uq, v, w, eta, mu, t);
      break;
    case 1:
      f = Source(x, uq, v, w, eta, mu, t);
      break;
    case 2:
      f = Tdfunc(x, uq, v, w, eta, mu, t);
      break;
    case 3:
      f = Fbou(x, uq, v, w, uhat, n, tau, eta, mu, t);
      break;
    case 4:
      f = Ubou(x, uq, v, w, uhat, n, tau, eta, mu, t);
      break;
    case 5:
      f = FbouHdg(x, uq, v, w, uhat, n, tau, eta, mu, t);
      break;
    case 6:
      f = Initu(x, eta, mu);
      break;
    case 7:
      f = VisScalars(x, uq, v, w, eta, mu, t);
      break;
    case 8:
      f = VisVectors(x, uq, v, w, eta, mu, t);
      break;
    case 9:
      f = QoIvolume(x, uq, v, w, eta, mu, t);
      break;
    case 10:
      f = QoIboundary(x, uq, v, w, uhat, n, tau, eta, mu, t);
      break;
    case 11:
      f = Fint(x, uq, v, w, uhat, n, tau, eta, mu, t);
      break;
    case 12:
      f = Fext(x, uq, v, w, uhat, n, uext, tau, eta, mu, t);
      break;
    default:
      throw std::runtime_error("Invalid function call in evaluateSymbolicFunctions");
  }

  return f;
}

void SymbolicScalarsVectors::func2cse(vec_pair &replacements, vec_basic &reduced_exprs, const std::vector<Expression> &f) {

   vec_basic exprs;
   for (const auto &fi : f) {
       exprs.push_back(fi.get_basic());
   }
   cse(replacements, reduced_exprs, exprs);
}

void SymbolicScalarsVectors::funcjac2cse(vec_pair &replacements,
                                         vec_basic &reduced_exprs_f,
                                         std::vector<vec_basic> &reduced_exprs_J,
                                         const std::vector<Expression> &f,
                                         const std::vector<std::vector<Expression>>& inputs_J) {
    vec_basic exprs;

    // Track original sizes
    int n_f = static_cast<int>(f.size());
    std::vector<int> jacobian_sizes;

    // Add original function expressions
    for (const auto &fi : f) {
        exprs.push_back(fi.get_basic());
    }

    // Append Jacobians and record sizes for each input group
    for (const auto& input_vec : inputs_J) {
        int m = f.size();
        int n = input_vec.size();
        jacobian_sizes.push_back(m * n);
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                exprs.push_back(f[i].diff(input_vec[j]).get_basic());
            }
        }
    }

    // Apply CSE
    vec_basic reduced_exprs;
    cse(replacements, reduced_exprs, exprs);

    // Decompose: f
    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + n_f);

    // Decompose: Jacobian blocks
    int offset = n_f;
    for (const int sz : jacobian_sizes) {
        vec_basic block(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_J.push_back(block);
        offset += sz;
    }
}

void SymbolicScalarsVectors::funcjachess2cse(vec_pair &replacements, vec_basic &reduced_exprs_f,
     std::vector<vec_basic> &reduced_exprs_J, std::vector<vec_basic> &reduced_exprs_H,
     const std::vector<Expression> &f, const std::vector<std::vector<Expression>>& inputs_J,
     const std::vector<std::vector<Expression>>& inputs_H) {

    vec_basic exprs;
    int count_f = f.size();
    int count_J = 0;
    int count_H = 0;

    // Add original function expressions
    for (const auto &fi : f) {
        exprs.push_back(fi.get_basic());
    }

    // Compute and append all Jacobian entries for each input group
    std::vector<int> J_sizes;
    for (const auto& input_vec : inputs_J) {
        int m = f.size();
        int n = input_vec.size();
        J_sizes.push_back(m * n);
        count_J += m * n;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                exprs.push_back(f[i].diff(input_vec[j]).get_basic());
            }
        }
    }

    // Compute and append all Hessian entries for each input group
    std::vector<int> H_sizes;
    for (const auto& input_vec : inputs_H) {
        int m = f.size();
        int n = input_vec.size();
        H_sizes.push_back(m * n * n);
        count_H += m * n * n;
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    exprs.push_back(f[i].diff(input_vec[j]).diff(input_vec[k]).get_basic());
                }
            }
        }
    }

    // Apply Common Subexpression Elimination
    vec_basic reduced_exprs;
    cse(replacements, reduced_exprs, exprs);

    // Decompose reduced_exprs into f, J, and H
    reduced_exprs_f.assign(reduced_exprs.begin(), reduced_exprs.begin() + count_f);

    int offset = count_f;
    for (int sz : J_sizes) {
        vec_basic Jblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_J.push_back(Jblock);
        offset += sz;
    }

    for (int sz : H_sizes) {
        vec_basic Hblock(reduced_exprs.begin() + offset, reduced_exprs.begin() + offset + sz);
        reduced_exprs_H.push_back(Hblock);
        offset += sz;
    }
}

void SymbolicScalarsVectors::func2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {
    std::ios_base::openmode mode = std::ios::out;
    if (append)
        mode |= std::ios::app;
    else
        mode |= std::ios::trunc;

    std::ofstream cppfile(filename + std::string(".cpp"), mode);
    cppfile << "void " <<funcname << "(dstype* f, ";
    cppfile << funcjacdecls[functionid] << "\n";
    cppfile << "{\n\n";

   if (f.size() > 0) {
       vec_pair replacements;
       vec_basic reduced_exprs;
       func2cse(replacements, reduced_exprs, f);

       // Determine variable usage
       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
       for (const auto &expr : f) {
           auto symbols = free_symbols(*expr.get_basic());
           used.insert(symbols.begin(), symbols.end());
       }

       auto depends_on = [&](const Expression &sym) {
           return used.count(sym.get_basic()) > 0;
       };

       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
       C99CodePrinter cpp;
       cppfile << "  Kokkos::parallel_for(\""<< funcnames[functionid] <<"\", N, KOKKOS_LAMBDA(const size_t i) {\n";
       
       // Emit symbolic variable loads
       for (const auto &[name, vec] : inputs) {
           for (size_t j = 0; j < vec.size(); ++j) {
               if (depends_on(vec[j])) { 
                 if (std::find(batch.begin(), batch.end(), name) != batch.end())
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
                 else 
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "];\n";
               }
           }
       }
       
       cppfile << "\n";
       
       // Emit intermediate CSE substitutions
       for (size_t n = 0; n < replacements.size(); ++n) {
           std::string var_name = cpp.apply(*replacements[n].first);
           std::string rhs = cpp.apply(*replacements[n].second);
           cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
       }
       cppfile << "\n";
       
       for (size_t n = 0; n < f.size(); ++n) {
           cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs[n]) << ";\n";
       }
       
       cppfile << "  });\n";
   }
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::funcjac2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {
    std::ios_base::openmode mode = std::ios::out;
    if (append)
        mode |= std::ios::app;
    else
        mode |= std::ios::trunc;

    std::ofstream cppfile(filename + std::string(".cpp"), mode);
    cppfile << "void " << funcname << "(dstype* f, ";
    int nJ = jacobianInputs[functionid].size();
    for (int k = 0; k < nJ; ++k)
        cppfile << "dstype* J" << (k+1) << ", ";
    cppfile << funcjacdecls[functionid] << "\n";
    cppfile << "{\n\n";

   if (f.size() > 0) {
       vec_pair replacements;
       vec_basic reduced_exprs_f;
       std::vector<vec_basic> reduced_exprs_J;
       funcjac2cse(replacements, reduced_exprs_f, reduced_exprs_J, f, jacobianInputs[functionid]);

       // Determine variable usage
       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
       for (const auto &expr : f) {
           auto symbols = free_symbols(*expr.get_basic());
           used.insert(symbols.begin(), symbols.end());
       }

       auto depends_on = [&](const Expression &sym) {
           return used.count(sym.get_basic()) > 0;
       };

       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
       C99CodePrinter cpp;
       cppfile << "  Kokkos::parallel_for(\""<< funcnames[functionid] <<"\", N, KOKKOS_LAMBDA(const size_t i) {\n";
       for (const auto &[name, vec] : inputs) {
           for (size_t j = 0; j < vec.size(); ++j) {
               if (depends_on(vec[j])) { 
                 if (std::find(batch.begin(), batch.end(), name) != batch.end())
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
                 else 
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "];\n";
               }
           }
       }

       for (size_t n = 0; n < replacements.size(); ++n) {
           std::string var_name = cpp.apply(*replacements[n].first);
           std::string rhs = cpp.apply(*replacements[n].second);
           cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
       }
       cppfile << "\n";
       
       for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {
           cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs_f[n]) << ";\n";
       }

       for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {
           for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {
               cppfile << "    J" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_J[k][j]) << ";\n";
           }
       }
       cppfile << "  });\n";
   }
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::funcjachess2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append) {
    vec_pair replacements;
    vec_basic reduced_exprs_f;
    std::vector<vec_basic> reduced_exprs_J, reduced_exprs_H;
    funcjachess2cse(replacements, reduced_exprs_f, reduced_exprs_J, reduced_exprs_H, f, jacobianInputs[functionid], hessianInputs[functionid]);

    // Determine variable usage
    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
    for (const auto &expr : f) {
        auto symbols = free_symbols(*expr.get_basic());
        used.insert(symbols.begin(), symbols.end());
    }

    auto depends_on = [&](const Expression &sym) {
        return used.count(sym.get_basic()) > 0;
    };

    std::ios_base::openmode mode = std::ios::out;
    if (append)
        mode |= std::ios::app;
    else
        mode |= std::ios::trunc;

    // Generate function prototype header based on functionid
    std::ofstream hfile(filename + std::string(".h"), mode);
    hfile << "void " << funcname << "jachess" << "(dstype* f, ";
    int nJ = reduced_exprs_J.size();
    int nH = reduced_exprs_H.size();
    for (int k = 0; k < nJ; ++k) hfile << "dstype* J" << (k+1) << ", ";
    for (int k = 0; k < nH; ++k) hfile << "dstype* H" << (k+1) << ", ";
    hfile << funcjacdecls[functionid] << ";\n";
    hfile.close();

    std::ofstream cppfile(filename + std::string(".cpp"), mode);
    cppfile << "void " << funcname << "jachess" << "(dstype* f, ";
    for (int k = 0; k < nJ; ++k) cppfile << "dstype* J" << (k+1) << ", ";
    for (int k = 0; k < nH; ++k) cppfile << "dstype* H" << (k+1) << ", ";
    cppfile << funcjacdecls[functionid] << "\n";
    cppfile << "{\n\n";

    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
    C99CodePrinter cpp;
       cppfile << "  Kokkos::parallel_for(\""<< funcnames[functionid] <<"\", N, KOKKOS_LAMBDA(const size_t i) {\n";
    for (const auto &[name, vec] : inputs) {
        for (size_t j = 0; j < vec.size(); ++j) {
            if (depends_on(vec[j])) { 
              if (std::find(batch.begin(), batch.end(), name) != batch.end())
                 cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*N+i];\n";
               else 
                 cppfile << "    dstype " << name << j << " = " << name << "[" << j << "];\n";
            }
        }
    }

    for (size_t n = 0; n < replacements.size(); ++n) {
        std::string var_name = cpp.apply(*replacements[n].first);
        std::string rhs = cpp.apply(*replacements[n].second);
        cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
    }
    cppfile << "\n";
    
    for (size_t n = 0; n < reduced_exprs_f.size(); ++n) {
        cppfile << "    f[" << n << " * N + i] = " << cpp.apply(*reduced_exprs_f[n]) << ";\n";
    }

    for (size_t k = 0; k < reduced_exprs_J.size(); ++k) {
        for (size_t j = 0; j < reduced_exprs_J[k].size(); ++j) {
            cppfile << "    J" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_J[k][j]) << ";\n";
        }
    }

    for (size_t k = 0; k < reduced_exprs_H.size(); ++k) {
        for (size_t j = 0; j < reduced_exprs_H[k].size(); ++j) {
            cppfile << "    H" << (k+1) << "[" << j << " * N + i] = " << cpp.apply(*reduced_exprs_H[k][j]) << ";\n";
        }
    }
       cppfile << "  });\n";
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::initfunc2cppfiles(const std::vector<Expression> &f, const std::string filename, const std::string funcname, const int functionid, bool append, int framework) {
    std::ios_base::openmode mode = std::ios::out;
    if (append)
        mode |= std::ios::app;
    else
        mode |= std::ios::trunc;

    std::ofstream cppfile(filename + std::string(".cpp"), mode);
    cppfile << "void " <<funcname << "(dstype* f, ";
    cppfile << "const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)" << "\n";
    cppfile << "{\n\n";

   if (f.size() > 0) {
       vec_pair replacements;
       vec_basic reduced_exprs;
       func2cse(replacements, reduced_exprs, f);

       // Determine variable usage
       std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
       for (const auto &expr : f) {
           auto symbols = free_symbols(*expr.get_basic());
           used.insert(symbols.begin(), symbols.end());
       }

       auto depends_on = [&](const Expression &sym) {
           return used.count(sym.get_basic()) > 0;
       };

       std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
       C99CodePrinter cpp;
       if (framework==0) 
       cppfile << "  for (int i = 0; i < N; ++i) {\n";
       else if (framework==1) 
       cppfile << "  Kokkos::parallel_for(\""<< funcnames[functionid] <<"\", N, KOKKOS_LAMBDA(const size_t i) {\n";
       else if (framework==2) 
       cppfile << "  parallel_for(\""<< funcnames[functionid] <<"\", N, GPU_LAMBDA(const size_t i) {\n";
       cppfile << "    int p = i%npe; \n";
       cppfile << "    int e = i/npe; \n";
       
       // Emit symbolic variable loads
       for (const auto &[name, vec] : inputs) {
           for (size_t j = 0; j < vec.size(); ++j) {
               if (depends_on(vec[j])) { 
                 if (std::find(batch.begin(), batch.end(), name) != batch.end())
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "*npe+p+npe*ncx*e];\n";
                 else 
                   cppfile << "    dstype " << name << j << " = " << name << "[" << j << "];\n";
               }
           }
       }
       
       cppfile << "\n";
       
       // Emit intermediate CSE substitutions
       for (size_t n = 0; n < replacements.size(); ++n) {
           std::string var_name = cpp.apply(*replacements[n].first);
           std::string rhs = cpp.apply(*replacements[n].second);
           cppfile << "    dstype " << var_name << " = " << rhs << ";\n";
       }
       cppfile << "\n";
       
       for (size_t n = 0; n < f.size(); ++n) {
           cppfile << "    f[p+npe*" << n << " +npe*nce*e] = " << cpp.apply(*reduced_exprs[n]) << ";\n";
       }
       
       if (framework==0) 
       cppfile << "  }\n";
       else if (framework==1) 
       cppfile << "  });\n";
       else if (framework==2) 
       cppfile << "  });\n";
   }
    cppfile << "}\n\n";
    cppfile.close();
}

void SymbolicScalarsVectors::appendUbouFbou(const std::string& filename, const std::string& funcname, int nbc) {
    std::ostringstream tmp;

    tmp << "void " << funcname << "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\n";
        tmp << "           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
        tmp << "           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\n";
        tmp << "           const int ncx, const int nco, const int ncw) {\n";

    for (int k = 1; k <= nbc; ++k) {
        if (k == 1)
            tmp << "    if (ib == 1 )\n";
        else
            tmp << "    else if (ib == " << k << " )\n";
        tmp << "        " << funcname << k << "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,\n";
        tmp << "                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);\n";
    }

    tmp << "}\n";

    std::ofstream cppfile(filename + ".cpp", std::ios::out | std::ios::app);
    cppfile << tmp.str();
    cppfile.close();
}
void SymbolicScalarsVectors::appendFbouHdg(const std::string& filename, const std::string& funcname, int nbc) {
    std::ostringstream tmp;

    tmp << "void " << funcname << "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\n";
        tmp << "           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
        tmp << "           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\n";
        tmp << "           const int ncx, const int nco, const int ncw) {\n";

    for (int k = 1; k <= nbc; ++k) {
        if (k == 1)
            tmp << "    if (ib == 1 )\n";
        else
            tmp << "    else if (ib == " << k << " )\n";
        tmp << "        " << funcname << k << "(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,\n";
        tmp << "                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);\n";
    }

    tmp << "}\n";

    std::ofstream cppfile(filename + ".cpp", std::ios::out | std::ios::app);
    cppfile << tmp.str();
    cppfile.close();
}
void SymbolicScalarsVectors::appendFextonly(const std::string& filename, const std::string& funcname, int nbc) {
    std::ostringstream tmp;

    tmp << "void " << funcname << "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\n";
        tmp << "           const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
        tmp << "           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\n";
        tmp << "           const int ncx, const int nco, const int ncw) {\n";

    for (int k = 1; k <= nbc; ++k) {
        if (k == 1)
            tmp << "    if (ib == 1 )\n";
        else
            tmp << "    else if (ib == " << k << " )\n";
        tmp << "        " << funcname << k << "(f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, modelnumber,\n";
        tmp << "                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd, ncx);\n";
    }

    tmp << "}\n";

    std::ofstream cppfile(filename + ".cpp", std::ios::out | std::ios::app);
    cppfile << tmp.str();
    cppfile.close();
}
void SymbolicScalarsVectors::appendFext(const std::string& filename, const std::string& funcname, int nbc) {
    std::ostringstream tmp;

    tmp << "void " << funcname << "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,\n";
        tmp << "           const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,\n";
        tmp << "           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,\n";
        tmp << "           const int ncx, const int nco, const int ncw) {\n";

    for (int k = 1; k <= nbc; ++k) {
        if (k == 1)
            tmp << "    if (ib == 1 )\n";
        else
            tmp << "    else if (ib == " << k << " )\n";
        tmp << "        " << funcname << k << "(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, modelnumber,\n";
        tmp << "                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd, ncx);\n";
    }

    tmp << "}\n";

    std::ofstream cppfile(filename + ".cpp", std::ios::out | std::ios::app);
    cppfile << tmp.str();
    cppfile.close();
}
void SymbolicScalarsVectors::emit_pointwise_value(
    std::ostream& os, const std::string& method_name,
    const std::string& cpp_signature,
    const std::vector<Expression>& f, int functionid)
{
    os << "    KOKKOS_INLINE_FUNCTION static\n";
    os << "    void " << method_name << "(" << cpp_signature << ") {\n";

    if (f.size() == 0) {
        os << "        // empty body — defaulted via ModelDefaults\n";
        os << "    }\n\n";
        return;
    }

    vec_pair replacements;
    vec_basic reduced_exprs;
    func2cse(replacements, reduced_exprs, f);

    std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
    for (const auto& expr : f) {
        auto symbols = free_symbols(*expr.get_basic());
        used.insert(symbols.begin(), symbols.end());
    }
    auto depends_on = [&](const Expression& sym) {
        return used.count(sym.get_basic()) > 0;
    };

    auto rename_input = [](const std::string& name) -> std::string {
        if (name == "eta")  return "uinf";
        if (name == "uhat") return "uh";
        return name;
    };

    std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
    C99CodePrinter cpp;
    for (const auto& [name, vec] : inputs) {
        for (size_t j = 0; j < vec.size(); ++j) {
            if (depends_on(vec[j])) {
                os << "        const double " << name << j
                   << " = " << rename_input(name) << "[" << j << "];\n";
            }
        }
    }

    if (!replacements.empty()) os << "\n";
    for (size_t n = 0; n < replacements.size(); ++n) {
        std::string var_name = cpp.apply(*replacements[n].first);
        std::string rhs      = cpp.apply(*replacements[n].second);
        os << "        const double " << var_name << " = " << rhs << ";\n";
    }
    os << "\n";

    for (size_t n = 0; n < f.size(); ++n) {
        os << "        f[" << n << "] = " << cpp.apply(*reduced_exprs[n]) << ";\n";
    }

    os << "    }\n\n";
}

void SymbolicScalarsVectors::emit_pointwise_value_per_ib(
    std::ostream& os, const std::string& method_name,
    const std::string& cpp_signature,
    const std::vector<Expression>& f, int functionid, int szuhat)
{
    os << "    KOKKOS_INLINE_FUNCTION static\n";
    os << "    void " << method_name << "(" << cpp_signature << ") {\n";

    int nbc = (szuhat > 0) ? (int)f.size() / szuhat : 0;
    if (nbc == 0) {
        os << "    }\n\n";
        return;
    }

    auto rename_input = [](const std::string& name) -> std::string {
        if (name == "eta")  return "uinf";
        if (name == "uhat") return "uh";
        return name;
    };

    C99CodePrinter cpp;
    for (int n = 0; n < nbc; ++n) {
        std::vector<Expression> g(szuhat);
        for (int m = 0; m < szuhat; ++m) g[m] = f[m + n * szuhat];

        os << "        " << ((n == 0) ? "if" : "else if")
           << " (ib == " << (n + 1) << ") {\n";

        vec_pair replacements;
        vec_basic reduced_exprs;
        func2cse(replacements, reduced_exprs, g);

        std::unordered_set<RCP<const Basic>, SymEngine::RCPBasicHash, SymEngine::RCPBasicKeyEq> used;
        for (const auto& expr : g) {
            auto symbols = free_symbols(*expr.get_basic());
            used.insert(symbols.begin(), symbols.end());
        }
        auto depends_on = [&](const Expression& sym) {
            return used.count(sym.get_basic()) > 0;
        };

        std::vector<std::pair<std::string, std::vector<Expression>>> inputs = inputvectors[functionid];
        for (const auto& [name, vec] : inputs) {
            for (size_t j = 0; j < vec.size(); ++j) {
                if (depends_on(vec[j])) {
                    os << "            const double " << name << j
                       << " = " << rename_input(name) << "[" << j << "];\n";
                }
            }
        }
        if (!replacements.empty()) os << "\n";
        for (size_t k = 0; k < replacements.size(); ++k) {
            std::string var_name = cpp.apply(*replacements[k].first);
            std::string rhs      = cpp.apply(*replacements[k].second);
            os << "            const double " << var_name << " = " << rhs << ";\n";
        }
        os << "\n";
        for (size_t k = 0; k < g.size(); ++k) {
            os << "            f[" << k << "] = " << cpp.apply(*reduced_exprs[k]) << ";\n";
        }
        os << "        }\n";
    }

    os << "    }\n\n";
}

void SymbolicScalarsVectors::generateModelHeader(const std::string& filename) {
    std::ofstream hfile(filename, std::ios::out | std::ios::trunc);

    hfile << "// Auto-generated by text2code. Do not edit by hand.\n";
    hfile << "// Regenerate with `text2code <pdeapp.txt>`.\n";
    hfile << "//\n";
    hfile << "// This header is consumed by `<exasim/model.hpp>`'s templated FEM\n";
    hfile << "// internals. The struct below satisfies the Model contract; the\n";
    hfile << "// inherited `ModelDefaults<GeneratedModel>` supplies zero-fill\n";
    hfile << "// defaults for any optional method this PDE doesn't define.\n";
    hfile << "#pragma once\n\n";
    hfile << "#include <Kokkos_Core.hpp>\n";
    hfile << "#include <exasim/model.hpp>\n\n";
    hfile << "struct GeneratedModel : exasim::ModelDefaults<GeneratedModel> {\n";
    hfile << "    static constexpr int nd     = 2;\n";
    hfile << "    static constexpr int ncu    = 1;\n";
    hfile << "    static constexpr int ncw    = 0;\n";
    hfile << "    static constexpr int nco    = 0;\n";
    hfile << "    static constexpr int nparam = 1;\n";
    hfile << "    static constexpr auto disc  = exasim::Discretization::HDG;\n";
    hfile << "    static constexpr int Nq = ncu * (1 + nd);\n\n";
    // ----- Volume value methods -----
    static const std::vector<std::tuple<std::string, std::string, std::string>>
        volume_methods = {
        {"Flux",       "flux",        "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
        {"Source",     "source",      "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
        {"Tdfunc",     "tdfunc",      "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
        {"VisScalars", "vis_scalars", "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
        {"VisVectors", "vis_vectors", "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
        {"QoIvolume",  "qoi_volume",  "double f[], const double x[], const double uq[], const double w[], const double mu[], const double uinf[], double t"},
    };
    for (const auto& [funcname, method_name, sig] : volume_methods) {
        auto it = std::find(funcnames.begin(), funcnames.end(), funcname);
        if (it == funcnames.end()) continue;
        int idx = it - funcnames.begin();
        if (!outputfunctions[idx]) continue;
        std::vector<Expression> f = evaluateSymbolicFunctions(idx);
        emit_pointwise_value(hfile, method_name, sig, f, idx);
    }

    // ----- Initial condition: initu(double ui[], const double x[], const double uinf[], const double mu[]) -----
    {
        auto it = std::find(funcnames.begin(), funcnames.end(), std::string("Initu"));
        if (it != funcnames.end()) {
            int idx = it - funcnames.begin();
            if (outputfunctions[idx]) {
                std::vector<Expression> f = evaluateSymbolicFunctions(idx);
                emit_pointwise_value(hfile, "initu",
                    "double f[], const double x[], const double uinf[], const double mu[]",
                    f, idx);
            }
        }
    }

    // ----- Boundary methods (per-ib dispatch) -----
    static const std::vector<std::tuple<std::string, std::string>>
        boundary_methods = {
        {"Fbou",        "fbou"},
        {"Ubou",        "ubou"},
        {"FbouHdg",     "fbou_hdg"},
        {"QoIboundary", "qoi_boundary"},
    };
    const std::string boundary_sig =
        "double f[], int ib, const double x[], const double uq[], const double w[],"
        " const double uh[], const double n[], const double tau[],"
        " const double mu[], const double uinf[], double t";
    for (const auto& [funcname, method_name] : boundary_methods) {
        auto it = std::find(funcnames.begin(), funcnames.end(), funcname);
        if (it == funcnames.end()) continue;
        int idx = it - funcnames.begin();
        if (!outputfunctions[idx]) continue;
        std::vector<Expression> f = evaluateSymbolicFunctions(idx);
        emit_pointwise_value_per_ib(hfile, method_name, boundary_sig, f, idx, szuhat);
    }

    // ----- HDG Jacobians (column-major: j outer, i inner) -----
    auto diff_to_exprs = [&](const std::vector<Expression>& f,
                             const std::vector<Expression>& input)
        -> std::vector<Expression> {
        std::vector<Expression> result;
        result.reserve(f.size() * input.size());
        for (size_t j = 0; j < input.size(); ++j) {
            for (size_t i = 0; i < f.size(); ++i) {
                result.emplace_back(f[i].diff(input[j]).get_basic());
            }
        }
        return result;
    };

    static const std::vector<std::tuple<std::string, std::string, std::vector<std::string>>>
        volume_jac_methods = {
        {"Flux",   "flux",   {"flux_jac_uq",   "flux_jac_w"}},
        {"Source", "source", {"source_jac_uq", "source_jac_w"}},
    };
    const std::string volume_sig =
        "double f[], const double x[], const double uq[], const double w[],"
        " const double mu[], const double uinf[], double t";
    for (const auto& [funcname, value_name, jac_names] : volume_jac_methods) {
        auto it = std::find(funcnames.begin(), funcnames.end(), funcname);
        if (it == funcnames.end()) continue;
        int idx = it - funcnames.begin();
        if (!outputfunctions[idx]) continue;
        std::vector<Expression> f = evaluateSymbolicFunctions(idx);
        // Skip Jacobians whose input vector is empty (e.g., w when ncw=0).
        const auto& jac_inputs = jacobianInputs[idx];
        for (size_t k = 0; k < jac_inputs.size() && k < jac_names.size(); ++k) {
            if (jac_inputs[k].empty()) continue;
            std::vector<Expression> jac = diff_to_exprs(f, jac_inputs[k]);
            emit_pointwise_value(hfile, jac_names[k], volume_sig, jac, idx);
        }
    }

    // FbouHdg Jacobians: per-ib dispatch, three Jacs (uq, w, uh).
    {
        auto it = std::find(funcnames.begin(), funcnames.end(), std::string("FbouHdg"));
        if (it != funcnames.end()) {
            int idx = it - funcnames.begin();
            if (outputfunctions[idx]) {
                std::vector<Expression> f = evaluateSymbolicFunctions(idx);
                int nbc = (szuhat > 0) ? (int)f.size() / szuhat : 0;
                const auto& jac_inputs = jacobianInputs[idx];
                static const std::vector<std::string> fbou_jac_names = {
                    "fbou_hdg_jac_uq", "fbou_hdg_jac_w", "fbou_hdg_jac_uh"};
                for (size_t k = 0; k < jac_inputs.size() && k < fbou_jac_names.size(); ++k) {
                    if (jac_inputs[k].empty()) continue;
                    // Build Jacobian per ib, concat into one vector with the
                    // same szuhat*njac slicing convention so the per-ib emitter
                    // can split it back out. Here we widen szuhat to the Jacobian
                    // block size for that input.
                    int jblock = szuhat * (int)jac_inputs[k].size();
                    std::vector<Expression> jac_all;
                    jac_all.reserve(jblock * nbc);
                    for (int n = 0; n < nbc; ++n) {
                        std::vector<Expression> g(szuhat);
                        for (int m = 0; m < szuhat; ++m) g[m] = f[m + n * szuhat];
                        std::vector<Expression> j_n = diff_to_exprs(g, jac_inputs[k]);
                        for (auto& e : j_n) jac_all.push_back(e);
                    }
                    emit_pointwise_value_per_ib(hfile, fbou_jac_names[k],
                        boundary_sig, jac_all, idx, jblock);
                }
            }
        }
    }

    hfile << "};\n";
    hfile.close();
}

