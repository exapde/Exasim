// #include <string>
// #include <vector>
// #include <iostream>
// #include <cmath>

#ifndef __COMPARESTRUCTS
#define __COMPARESTRUCTS

// ---- Helper comparison functions ----

static void printElemRecv(const std::vector<std::array<int,3>>& elemrecv,
                          const std::string& name = "elemrecv")
{
    std::cout << name << " (size = " << elemrecv.size() << "):\n";

    for (std::size_t i = 0; i < elemrecv.size(); ++i) {
        const auto& e = elemrecv[i];
        std::cout << e[0] << "  " << e[1] << "  " << e[2] << " \n";
    }
}

static bool compareString(const std::string& a,
                          const std::string& b,
                          const std::string& name,
                          bool verbose)
{
    if (a != b) {
        if (verbose) {
            std::cout << "Mismatch in " << name << ": \""
                      << a << "\" vs \"" << b << "\"\n";
        }
        return false;
    }
    return true;
}

static bool compareVecInt(const std::vector<int>& a,
                          const std::vector<int>& b,
                          const std::string& name,
                          bool verbose)
{
    if (a.size() != b.size()) {
        if (verbose) {
            std::cout << "Mismatch in " << name << " size: "
                      << a.size() << " vs " << b.size() << "\n";
        }
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            if (verbose) {
                std::cout << "Mismatch in " << name << "[" << i << "]: "
                          << a[i] << " vs " << b[i] << "\n";
            }
            return false;
        }
    }
    return true;
}

static bool compareVecIntSorted(const std::vector<int>& a,
                                const std::vector<int>& b,
                                const std::string& name,
                                bool verbose)
{
    // First check size; no need to sort if sizes differ
    if (a.size() != b.size()) {
        if (verbose) {
            std::cout << "Mismatch in " << name << " size: "
                      << a.size() << " vs " << b.size() << "\n";
        }
        return false;
    }

    // Make sorted copies
    std::vector<int> sa = a;
    std::vector<int> sb = b;
    std::sort(sa.begin(), sa.end());
    std::sort(sb.begin(), sb.end());

    // Compare sorted arrays
    for (std::size_t i = 0; i < sa.size(); ++i) {
        if (sa[i] != sb[i]) {
            if (verbose) {
                std::cout << "Mismatch in sorted " << name << "[" << i << "]: "
                          << sa[i] << " vs " << sb[i] << "\n";
            }
            return false;
        }
    }

    return true;
}

static bool compareVecDouble(const std::vector<double>& a,
                             const std::vector<double>& b,
                             const std::string& name,
                             bool verbose,
                             double tol)
{
    if (a.size() != b.size()) {
        if (verbose) {
            std::cout << "Mismatch in " << name << " size: "
                      << a.size() << " vs " << b.size() << "\n";
        }
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        double diff = std::abs(a[i] - b[i]);
        if (tol == 0.0) {
            if (a[i] != b[i]) {
                if (verbose) {
                    std::cout << "Mismatch in " << name << "[" << i << "]: "
                              << a[i] << " vs " << b[i] << "\n";
                }
                return false;
            }
        } else {
            if (diff > tol) {
                if (verbose) {
                    std::cout << "Mismatch in " << name << "[" << i << "]: "
                              << a[i] << " vs " << b[i]
                              << " (|diff| = " << diff
                              << ", tol = " << tol << ")\n";
                }
                return false;
            }
        }
    }
    return true;
}

static bool compareDouble(double a, double b,
                          const std::string& name,
                          bool verbose,
                          double tol)
{
    double diff = std::abs(a - b);
    if (tol == 0.0) {
        if (a != b) {
            if (verbose) {
                std::cout << "Mismatch in " << name << ": "
                          << a << " vs " << b << "\n";
            }
            return false;
        }
    } else {
        if (diff > tol) {
            if (verbose) {
                std::cout << "Mismatch in " << name << ": "
                          << a << " vs " << b
                          << " (|diff| = " << diff
                          << ", tol = " << tol << ")\n";
            }
            return false;
        }
    }
    return true;
}

static bool compareArray3(const std::vector<std::array<int,3>>& a,
                          const std::vector<std::array<int,3>>& b,
                          const std::string& name,
                          bool verbose)
{
    if (a.size() != b.size()) {
        if (verbose)
            std::cout << "Mismatch in " << name << ": size "
                      << a.size() << " vs " << b.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            if (verbose) {
                std::cout << "Mismatch in " << name << "[" << i << "]: "
                          << "[" << a[i][0] << "," << a[i][1] << "," << a[i][2] << "] vs "
                          << "[" << b[i][0] << "," << b[i][1] << "," << b[i][2] << "]\n";
            }
            return false;
        }
    }
    return true;
}

// ---- Main comparison function ----

void comparePDE(const PDE& p1,
                const PDE& p2,
                bool verbose = false,
                double tol = 0.0)
{
    bool all_equal = true;  // internal flag, not returned

    // Strings
    if (!compareString(p1.discretization, p2.discretization, "discretization", verbose))
        all_equal = false;
    if (!compareString(p1.exasimpath,     p2.exasimpath,     "exasimpath",     verbose))
        all_equal = false;
    if (!compareString(p1.datapath,       p2.datapath,       "datapath",       verbose))
        all_equal = false;
    if (!compareString(p1.datainpath,     p2.datainpath,     "datainpath",     verbose))
        all_equal = false;
    if (!compareString(p1.dataoutpath,    p2.dataoutpath,    "dataoutpath",    verbose))
        all_equal = false;
    if (!compareString(p1.platform,       p2.platform,       "platform",       verbose))
        all_equal = false;
    if (!compareString(p1.model,          p2.model,          "model",          verbose))
        all_equal = false;
    if (!compareString(p1.pdeappfile,     p2.pdeappfile,     "pdeappfile",     verbose))
        all_equal = false;
    if (!compareString(p1.modelfile,      p2.modelfile,      "modelfile",      verbose))
        all_equal = false;
    if (!compareString(p1.meshfile,       p2.meshfile,       "meshfile",       verbose))
        all_equal = false;
    if (!compareString(p1.xdgfile,        p2.xdgfile,        "xdgfile",        verbose))
        all_equal = false;
    if (!compareString(p1.udgfile,        p2.udgfile,        "udgfile",        verbose))
        all_equal = false;
    if (!compareString(p1.vdgfile,        p2.vdgfile,        "vdgfile",        verbose))
        all_equal = false;
    if (!compareString(p1.wdgfile,        p2.wdgfile,        "wdgfile",        verbose))
        all_equal = false;
    if (!compareString(p1.uhatfile,       p2.uhatfile,       "uhatfile",       verbose))
        all_equal = false;
    if (!compareString(p1.partitionfile,  p2.partitionfile,  "partitionfile",  verbose))
        all_equal = false;

    // Ints – no early return; just mark mismatch
    #define CMP_INT(field) \
        do { \
            if (p1.field != p2.field) { \
                if (verbose) { \
                    std::cout << "Mismatch in " #field ": " \
                              << p1.field << " vs " << p2.field << "\n"; \
                } \
                all_equal = false; \
            } \
        } while (0)

    CMP_INT(gencode);
    CMP_INT(writemeshsol);
    CMP_INT(modelnumber);
    CMP_INT(mpiprocs);
    CMP_INT(nd);
    CMP_INT(nc);
    CMP_INT(ncu);
    CMP_INT(ncq);
    CMP_INT(ncp);
    CMP_INT(ncv);
    CMP_INT(nch);
    CMP_INT(ncx);
    CMP_INT(ncw);
    CMP_INT(nce);
    CMP_INT(np);
    CMP_INT(nve);
    CMP_INT(ne);
    CMP_INT(nsca);
    CMP_INT(nvec);
    CMP_INT(nten);
    CMP_INT(nsurf);
    CMP_INT(nvqoi);
    CMP_INT(neb);
    CMP_INT(nfb);
    CMP_INT(elemtype);
    CMP_INT(nodetype);
    CMP_INT(hybrid);
    CMP_INT(tdep);
    CMP_INT(wave);
    CMP_INT(linearproblem);
    CMP_INT(subproblem);
    CMP_INT(debugmode);
    CMP_INT(stgNmode);
    CMP_INT(porder);
    CMP_INT(pgauss);
    CMP_INT(temporalscheme);
    CMP_INT(torder);
    CMP_INT(nstage);
    CMP_INT(convStabMethod);
    CMP_INT(diffStabMethod);
    CMP_INT(rotatingFrame);
    CMP_INT(viscosityModel);
    CMP_INT(SGSmodel);
    CMP_INT(ALE);
    CMP_INT(AV);
    CMP_INT(AVdistfunction);
    CMP_INT(AVsmoothingIter);
    CMP_INT(frozenAVflag);
    CMP_INT(nonlinearsolver);
    CMP_INT(linearsolver);
    CMP_INT(NewtonIter);
    CMP_INT(GMRESiter);
    CMP_INT(GMRESrestart);
    CMP_INT(GMRESortho);
    CMP_INT(preconditioner);
    CMP_INT(precMatrixType);
    CMP_INT(ppdegree);
    CMP_INT(NLMatrixType);
    CMP_INT(runmode);
    CMP_INT(tdfunc);
    CMP_INT(sourcefunc);
    CMP_INT(matvecorder);
    CMP_INT(RBdim);
    CMP_INT(saveSolFreq);
    CMP_INT(saveSolOpt);
    CMP_INT(timestepOffset);
    CMP_INT(saveSolBouFreq);
    CMP_INT(ibs);
    CMP_INT(compudgavg);
    CMP_INT(extFhat);
    CMP_INT(extUhat);
    CMP_INT(extStab);
    CMP_INT(saveResNorm);
    CMP_INT(dae_steps);
    CMP_INT(coupledinterface);
    CMP_INT(coupledcondition);
    CMP_INT(coupledboundarycondition);

    #undef CMP_INT

    // Doubles (with tolerance)
    if (!compareDouble(p1.time,        p2.time,        "time",        verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.NLparam,     p2.NLparam,     "NLparam",     verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.NewtonTol,   p2.NewtonTol,   "NewtonTol",   verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.GMREStol,    p2.GMREStol,    "GMREStol",    verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.matvectol,   p2.matvectol,   "matvectol",   verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.dae_alpha,   p2.dae_alpha,   "dae_alpha",   verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.dae_beta,    p2.dae_beta,    "dae_beta",    verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.dae_gamma,   p2.dae_gamma,   "dae_gamma",   verbose, tol))
        all_equal = false;
    if (!compareDouble(p1.dae_epsilon, p2.dae_epsilon, "dae_epsilon", verbose, tol))
        all_equal = false;

    // Vectors<int>
    if (!compareVecInt(p1.interfaceFluxmap, p2.interfaceFluxmap, "interfaceFluxmap", verbose))
        all_equal = false;

    // Vectors<double>
    if (!compareVecDouble(p1.dae_dt,       p2.dae_dt,       "dae_dt",       verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.dt,           p2.dt,           "dt",           verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.tau,          p2.tau,          "tau",          verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.flag,         p2.flag,         "flag",         verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.problem,      p2.problem,      "problem",      verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.solversparam, p2.solversparam, "solversparam", verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.factor,       p2.factor,       "factor",       verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.physicsparam, p2.physicsparam, "physicsparam", verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.externalparam,p2.externalparam,"externalparam",verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.vindx,        p2.vindx,        "vindx",        verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.avparam1,     p2.avparam1,     "avparam1",     verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.avparam2,     p2.avparam2,     "avparam2",     verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.stgib,        p2.stgib,        "stgib",        verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.stgdata,      p2.stgdata,      "stgdata",      verbose, tol))
        all_equal = false;
    if (!compareVecDouble(p1.stgparam,     p2.stgparam,     "stgparam",     verbose, tol))
        all_equal = false;

    // Optional summary
    if (verbose) {
        if (all_equal)
            std::cout << "PDE instances are identical.\n";
        else
            std::cout << "PDE instances differ.\n";
    }

    // no return value
}

void compareMaster(const Master& m1,
                   const Master& m2,
                   bool verbose = false,
                   double tol = 0.0)
{
    bool all_equal = true;  // internal flag only for summary

    // ---- vector<double> fields ----
    #define CMP_VECD(field) \
        do { \
            if (!compareVecDouble(m1.field, m2.field, #field, verbose, tol)) \
                all_equal = false; \
        } while (0)

    CMP_VECD(xpe);
    CMP_VECD(gpe);
    CMP_VECD(gwe);
    CMP_VECD(xpf);
    CMP_VECD(gpf);
    CMP_VECD(gwf);

    CMP_VECD(shapeg);
    CMP_VECD(shapegt);
    CMP_VECD(shapfg);
    CMP_VECD(shapfgt);

    CMP_VECD(shapent);
    CMP_VECD(shapfnt);
    CMP_VECD(shapegw);
    CMP_VECD(shapfgw);

    CMP_VECD(shapen);
    CMP_VECD(shapfn);

    CMP_VECD(xp1d);
    CMP_VECD(gp1d);
    CMP_VECD(gw1d);

    CMP_VECD(shap1dg);
    CMP_VECD(shap1dgt);
    CMP_VECD(shap1dn);
    CMP_VECD(shap1dnt);
    CMP_VECD(shap1dgw);

    CMP_VECD(phielem);
    CMP_VECD(phiface);

    #undef CMP_VECD

    // ---- vector<int> fields ----
    #define CMP_VECI(field) \
        do { \
            if (!compareVecInt(m1.field, m2.field, #field, verbose)) \
                all_equal = false; \
        } while (0)

    CMP_VECI(telem);
    CMP_VECI(tface);
    CMP_VECI(perm);
    CMP_VECI(permind);

    #undef CMP_VECI

    // ---- scalar int fields ----
    #define CMP_INT(field) \
        do { \
            if (m1.field != m2.field) { \
                if (verbose) \
                    std::cout << "Mismatch in " #field ": " \
                              << m1.field << " vs " << m2.field << "\n"; \
                all_equal = false; \
            } \
        } while (0)

    CMP_INT(nd);
    CMP_INT(npe);
    CMP_INT(npf);
    CMP_INT(nge);
    CMP_INT(ngf);
    CMP_INT(porder);
    CMP_INT(pgauss);
    CMP_INT(nfe);
    CMP_INT(elemtype);
    CMP_INT(nodetype);
    CMP_INT(nve);
    CMP_INT(nvf);
    CMP_INT(np1d);
    CMP_INT(ng1d);
    CMP_INT(npermind);

    #undef CMP_INT

    // ---- optional summary ----
    if (verbose) {
        if (all_equal)
            std::cout << "Master instances are identical.\n";
        else
            std::cout << "Master instances differ.\n";
    }
}

void compareDMD(const DMD& d1, const DMD& d2, bool verbose = false)
{
    bool all_equal = true;

    // ---- vector<int> fields ----
    #define CMP_VECI(field) \
        do { \
            if (!compareVecInt(d1.field, d2.field, #field, verbose)) \
                all_equal = false; \
        } while (0)

    CMP_VECI(nbsd);
    CMP_VECI(elempart);
    CMP_VECI(elem2cpu);
    CMP_VECI(elemsendpts);
    CMP_VECI(elemrecvpts);
    CMP_VECI(elempartpts);  
    CMP_VECI(intepartpts);
    //CMP_VECI(nbinfo);

    #undef CMP_VECI

    // ---- vector<array<int,3>> fields ----
    #define CMP_ARR3(field) \
        do { \
            if (!compareArray3(d1.field, d2.field, #field, verbose)) \
                all_equal = false; \
        } while (0)

    CMP_ARR3(elemrecv);
    CMP_ARR3(elemsend);

    #undef CMP_ARR3

    // ---- scalar int field ----
    // if (d1.numneigh != d2.numneigh) {
    //     if (verbose) {
    //         std::cout << "Mismatch in numneigh: "
    //                   << d1.numneigh << " vs " << d2.numneigh << "\n";
    //     }
    //     all_equal = false;
    // }

    // ---- optional summary ----
    if (verbose) {
        if (all_equal)
            std::cout << "DMD instances are identical.\n";
        else
            std::cout << "DMD instances differ.\n";
    }
}

#endif

// static bool compareVector(const std::vector<int>& a,
//                           const std::vector<int>& b,
//                           const std::string& name,
//                           bool verbose)
// {
//     if (a.size() != b.size()) {
//         if (verbose)
//             std::cout << "Mismatch in " << name << ": size "
//                       << a.size() << " vs " << b.size() << "\n";
//         return false;
//     }
//     for (size_t i = 0; i < a.size(); ++i) {
//         if (a[i] != b[i]) {
//             if (verbose)
//                 std::cout << "Mismatch in " << name << "[" << i << "]: "
//                           << a[i] << " vs " << b[i] << "\n";
//             return false;
//         }
//     }
//     return true;
// }

// static bool compareVecDouble(const std::vector<double>& a,
//                              const std::vector<double>& b,
//                              const std::string& name,
//                              bool verbose,
//                              double tol = 0.0) // exact compare by default
// {
//     if (a.size() != b.size()) {
//         if (verbose)
//             std::cout << "Mismatch in " << name << " size: "
//                       << a.size() << " vs " << b.size() << "\n";
//         return false;
//     }
// 
//     for (size_t i = 0; i < a.size(); ++i) {
//         if (tol == 0.0) {
//             if (a[i] != b[i]) {
//                 if (verbose)
//                     std::cout << "Mismatch in " << name << "[" << i << "]: "
//                               << a[i] << " vs " << b[i] << "\n";
//                 return false;
//             }
//         } else {
//             if (std::abs(a[i] - b[i]) > tol) {
//                 if (verbose)
//                     std::cout << "Mismatch in " << name << "[" << i << "]: "
//                               << a[i] << " vs " << b[i]
//                               << " (tol = " << tol << ")\n";
//                 return false;
//             }
//         }
//     }
// 
//     return true;
// }
// 
// static bool compareVecInt(const std::vector<int>& a,
//                           const std::vector<int>& b,
//                           const std::string& name,
//                           bool verbose)
// {
//     if (a.size() != b.size()) {
//         if (verbose)
//             std::cout << "Mismatch in " << name << " size: "
//                       << a.size() << " vs " << b.size() << "\n";
//         return false;
//     }
// 
//     for (size_t i = 0; i < a.size(); ++i) {
//         if (a[i] != b[i]) {
//             if (verbose)
//                 std::cout << "Mismatch in " << name << "[" << i << "]: "
//                           << a[i] << " vs " << b[i] << "\n";
//             return false;
//         }
//     }
//     return true;
// }

// // Helper to compare vectors
// template<typename T>
// bool vec_equal(const std::vector<T>& a, const std::vector<T>& b) {
//     return a == b;
// }
// 
// // Helper to compare 2D vectors
// template<typename T>
// bool vec2_equal(const std::vector<std::vector<T>>& a,
//                 const std::vector<std::vector<T>>& b) {
//     if (a.size() != b.size()) return false;
//     for (size_t i = 0; i < a.size(); i++) {
//         if (a[i] != b[i]) return false;
//     }
//     return true;
// }
// 
// // Main comparison function
// bool compareDMD(const DMD& d1, const DMD& d2, bool verbose = false) {
//     if (!vec_equal(d1.nbsd, d2.nbsd)) {
//         if (verbose) std::cout << "Mismatch: nbsd\n";
//         return false;
//     }
//     if (!vec2_equal(d1.elemrecv, d2.elemrecv)) {
//         if (verbose) std::cout << "Mismatch: elemrecv\n";
//         return false;
//     }
//     if (!vec2_equal(d1.elemsend, d2.elemsend)) {
//         if (verbose) std::cout << "Mismatch: elemsend\n";
//         return false;
//     }
//     if (!vec_equal(d1.elempart, d2.elempart)) {
//         if (verbose) std::cout << "Mismatch: elempart\n";
//         return false;
//     }
//     if (!vec_equal(d1.elem2cpu, d2.elem2cpu)) {
//         if (verbose) std::cout << "Mismatch: elem2cpu\n";
//         return false;
//     }
//     if (!vec_equal(d1.elemsendpts, d2.elemsendpts)) {
//         if (verbose) std::cout << "Mismatch: elemsendpts\n";
//         return false;
//     }
//     if (!vec_equal(d1.elemrecvpts, d2.elemrecvpts)) {
//         if (verbose) std::cout << "Mismatch: elemrecvpts\n";
//         return false;
//     }
//     if (!vec_equal(d1.elempartpts, d2.elempartpts)) {
//         if (verbose) std::cout << "Mismatch: elempartpts\n";
//         return false;
//     }
//     if (!vec_equal(d1.intepartpts, d2.intepartpts)) {
//         if (verbose) std::cout << "Mismatch: intepartpts\n";
//         return false;
//     }
//     if (!vec_equal(d1.nbinfo, d2.nbinfo)) {
//         if (verbose) std::cout << "Mismatch: nbinfo\n";
//         return false;
//     }
//     if (d1.numneigh != d2.numneigh) {
//         if (verbose) std::cout << "Mismatch: numneigh\n";
//         return false;
//     }
// 
//     return true; // All fields match
// }