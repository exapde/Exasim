#include "SymbolicFunctions.cpp"
#include "SymbolicScalarsVectors.cpp"

int main() 
{
  SymbolicScalarsVectors ssv;

  for (int i=0; i<ssv.outputfunctions.size(); i++) {
    std::string funcname = ssv.funcnames[i];
    if (ssv.outputfunctions[i] == true) {
      std::vector<Expression> f = ssv.evaluateSymbolicFunctions(i);
      bool append = true;
      if (ssv.exasim == true) {
        std::string fname = std::string("Kokkos")  + funcname;
        std::string jname = std::string("Hdg") + funcname;
        if (funcname == "FbouHdg") {
          fname = std::string("HdgFbouonly");
          jname = std::string("HdgFbou");
        }
        if (funcname == "Sourcew") {
          fname = std::string("HdgSourcewonly");
          jname = std::string("HdgSourcew");
        }
        if ((funcname == "Ubou") || (funcname == "Fbou") || (funcname == "FbouHdg")) { 
          int szf = f.size();
          int szuhat = ssv.szuhat;
          int nbc = szf/szuhat;
          for (int n = 0; n < nbc; ++n) {
            std::vector<Expression> g(szuhat);
            for (int m = 0; m < szuhat; ++m) {
              g[m] = f[m + n * szuhat];
            }
            if (n==0) {
               ssv.func2cppfiles(g, ssv.modelpath + fname, fname + std::to_string(n+1), i, false);
               if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(g, ssv.modelpath + jname, jname + std::to_string(n+1), i, false);
            } else {
               ssv.func2cppfiles(g, ssv.modelpath + fname, fname + std::to_string(n+1), i, append);
               if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(g, ssv.modelpath + jname, jname + std::to_string(n+1), i, append);
            }
            if (n==nbc-1) {
               ssv.appendUbouFbou(ssv.modelpath + fname, fname, nbc);
               if (funcname == "FbouHdg") ssv.appendFbouHdg(ssv.modelpath + jname, jname, nbc);
            }
          }
        } else if ((funcname == "Initu") || (funcname == "Initq") || (funcname == "Inituq") || (funcname == "Initv") || (funcname == "Initw")) { 
          std::string kname = funcname;
          if (funcname == "Inituq") kname = "Initudg";
          if (funcname == "Initv") kname = "Initvdg";
          if (funcname == "Initw") kname = "Initwdg";
          ssv.initfunc2cppfiles(f, ssv.modelpath + "cpu" + kname, "cpu" + kname, i, false, 0);
          ssv.initfunc2cppfiles(f, ssv.modelpath + "Kokkos" + kname, "Kokkos" + kname, i, false, 1);
        } else {
          ssv.func2cppfiles(f, ssv.modelpath + fname, fname, i, false);
          if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, ssv.modelpath + jname, jname, i, false);
        }
      } else {
        ssv.func2cppfiles(f, funcname, funcname, i, false);
        if (ssv.jacobianInputs[i].size() > 0) ssv.funcjac2cppfiles(f, funcname, funcname, i, append);
        if (ssv.hessianInputs[i].size() > 0) ssv.funcjachess2cppfiles(f, funcname, funcname, i, append);
      }
    }
  }
}
