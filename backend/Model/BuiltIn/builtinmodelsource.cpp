#include <cstdio>
#include <cstdlib>

#include "model1/model.cpp"
#include "model2/model.cpp"
#include "model3/model.cpp"
#include "model4/model.cpp"

namespace m1 = exasim_model_1;
namespace m2 = exasim_model_2;
namespace m3 = exasim_model_3;
namespace m4 = exasim_model_4;

// ----------------------------- Kokkos -----------------------------

extern "C"
void builtinKokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uinf, const dstype* param,
                dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosFlux(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosFlux(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;      
      case 3: m3::KokkosFlux(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosFlux(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosFlux\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosAvfield(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, const dstype* uinf, const dstype* param,
                   dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                   int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosAvfield(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosAvfield(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;      
      case 3: m3::KokkosAvfield(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosAvfield(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosAvfield\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosEoS(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
               const dstype* wdg, const dstype* uinf, const dstype* param,
               dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
               int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosEoS(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosEoS(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 3: m3::KokkosEoS(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosEoS(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosEoS\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosEoSdu(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uinf, const dstype* param,
                 dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                 int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosEoSdu(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosEoSdu(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;     
      case 3: m3::KokkosEoSdu(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosEoSdu(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosEoSdu\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosEoSdw(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uinf, const dstype* param,
                 dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                 int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosEoSdw(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosEoSdw(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;    
      case 3: m3::KokkosEoSdw(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosEoSdw(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosEoSdw\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau,
                const dstype* uinf, const dstype* param, dstype time, int builtinmodelID,
                int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosFbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosFbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return; 
      case 3: m3::KokkosFbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosFbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosFbou\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosFhat(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,
                const dstype* odg1, const dstype* odg2, const dstype* wdg1, const dstype* wdg2,
                const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int builtinmodelID, int ng, int nc, int ncu,
                int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosFhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosFhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;     
      case 3: m3::KokkosFhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosFhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosFhat\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                   int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::KokkosInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;   
      case 3: m3::KokkosInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::KokkosInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosInitodg\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                 int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::KokkosInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;      
      case 3: m3::KokkosInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::KokkosInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosInitq\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                 int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::KokkosInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;     
      case 3: m3::KokkosInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::KokkosInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosInitu\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                   int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::KokkosInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;    
      case 3: m3::KokkosInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::KokkosInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosInitudg\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                   int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::KokkosInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 3: m3::KokkosInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::KokkosInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosInitwdg\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosMonitor(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, const dstype* uinf, const dstype* param,
                   dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                   int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosMonitor(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosMonitor(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;    
      case 3: m3::KokkosMonitor(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosMonitor(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosMonitor\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosOutput(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                  const dstype* wdg, const dstype* uinf, const dstype* param,
                  dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                  int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosOutput(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosOutput(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;    
      case 3: m3::KokkosOutput(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosOutput(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosOutput\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                  const dstype* wdg, const dstype* uinf, const dstype* param,
                  dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                  int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosSource(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosSource(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;  
      case 3: m3::KokkosSource(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosSource(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosSource\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosSourcew(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, const dstype* uinf, const dstype* param,
                   dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                   int ncx, int nco, int ncw, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosSourcew(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 2: m2::KokkosSourcew(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;      
      case 3: m3::KokkosSourcew(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      case 4: m4::KokkosSourcew(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosSourcew\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosStab(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,
                const dstype* odg1, const dstype* odg2, const dstype* wdg1, const dstype* wdg2,
                const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int builtinmodelID, int ng, int nc, int ncu,
                int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosStab(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosStab(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;      
      case 3: m3::KokkosStab(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosStab(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosStab\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosTdfunc(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                  const dstype* wdg, const dstype* uinf, const dstype* param,
                  dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                  int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosTdfunc(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosTdfunc(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;    
      case 3: m3::KokkosTdfunc(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosTdfunc(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosTdfunc\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosUbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau,
                const dstype* uinf, const dstype* param, dstype time, int builtinmodelID, int ib,
                int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosUbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosUbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;     
      case 3: m3::KokkosUbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosUbou(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosUbou\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosUhat(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,
                const dstype* odg1, const dstype* odg2, const dstype* wdg1, const dstype* wdg2,
                const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int builtinmodelID, int ng, int nc, int ncu,
                int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosUhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosUhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;    
      case 3: m3::KokkosUhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosUhat(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosUhat\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uinf, const dstype* param,
                      dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosVisScalars(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosVisScalars(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::KokkosVisScalars(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosVisScalars(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosVisScalars\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosVisVectors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uinf, const dstype* param,
                      dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosVisVectors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosVisVectors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::KokkosVisVectors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosVisVectors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosVisVectors\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosVisTensors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uinf, const dstype* param,
                      dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosVisTensors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosVisTensors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::KokkosVisTensors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosVisTensors(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosVisTensors\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosQoIvolume(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                     const dstype* wdg, const dstype* uinf, const dstype* param,
                     dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd,
                     int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosQoIvolume(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosQoIvolume(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::KokkosQoIvolume(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosQoIvolume(f, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosQoIvolume\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinKokkosQoIboundary(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                       const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau,
                       const dstype* uinf, const dstype* param, dstype time, int builtinmodelID,
                       int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::KokkosQoIboundary(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::KokkosQoIboundary(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;  
      case 3: m3::KokkosQoIboundary(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::KokkosQoIboundary(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in KokkosQoIboundary\n", builtinmodelID); std::abort();
    }
}

// ----------------------------- cpuInit* -----------------------------

extern "C"
void builtincpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::cpuInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::cpuInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;     
      case 3: m3::cpuInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::cpuInitodg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in cpuInitodg\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtincpuInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
              int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::cpuInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::cpuInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;     
      case 3: m3::cpuInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::cpuInitq(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in cpuInitq\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtincpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
              int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::cpuInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::cpuInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return; 
      case 3: m3::cpuInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::cpuInitu(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in cpuInitu\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtincpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::cpuInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::cpuInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;    
      case 3: m3::cpuInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::cpuInitudg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in cpuInitudg\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtincpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param,
                int builtinmodelID, int ng, int ncx, int nce, int npe, int ne)
{
    switch (builtinmodelID) {
      case 1: m1::cpuInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 2: m2::cpuInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;    
      case 3: m3::cpuInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      case 4: m4::cpuInitwdg(f, xdg, uinf, param, builtinmodelID, ng, ncx, nce, npe, ne); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in cpuInitwdg\n", builtinmodelID); std::abort();
    }
}

// ----------------------------- HDG -----------------------------

extern "C"
void builtinHdgEoS(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg,
            const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param,
            dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd, int ncx,
            int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgEoS(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgEoS(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;  
      case 3: m3::HdgEoS(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgEoS(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgEoS\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg,
             const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
             const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param,
             dstype time, int builtinmodelID, int ib, int ng, int nc, int ncu, int nd,
             int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFbou(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFbou(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;    
      case 3: m3::HdgFbou(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFbou(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFbou\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau,
                 const dstype* uinf, const dstype* param, dstype time, int builtinmodelID,
                 int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFbouonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFbouonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return; 
      case 3: m3::HdgFbouonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFbouonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFbouonly\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg,
             const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
             const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param,
             dstype time, int builtinmodelID, int ib, int ng, int nc, int ncu, int nd,
             int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFint(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFint(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::HdgFint(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFint(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFint\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFintonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau,
                 const dstype* uinf, const dstype* param, dstype time, int builtinmodelID,
                 int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFintonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFintonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;      
      case 3: m3::HdgFintonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFintonly(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFintonly\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFext(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg,
             const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
             const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf,
             const dstype* param, dstype time, int builtinmodelID, int ib, int ng, int nc,
             int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFext(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFext(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::HdgFext(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFext(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFext\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFextonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* uext,
                 const dstype* tau, const dstype* uinf, const dstype* param,
                 dstype time, int builtinmodelID, int ib, int ng, int nc, int ncu, int nd,
                 int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFextonly(f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFextonly(f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;     
      case 3: m3::HdgFextonly(f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFextonly(f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, builtinmodelID, ib, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFextonly\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg,
             const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param,
             dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd, int ncx,
             int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgFlux(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgFlux(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;      
      case 3: m3::HdgFlux(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgFlux(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgFlux\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg,
               const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param,
               dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd, int ncx,
               int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgSource(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgSource(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;     
      case 3: m3::HdgSource(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgSource(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgSource\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgSourcew(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param,
                dstype time, int builtinmodelID, int ng, int nc, int ncu, int nd, int ncx,
                int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgSourcew(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgSourcew(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;   
      case 3: m3::HdgSourcew(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgSourcew(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgSourcew\n", builtinmodelID); std::abort();
    }
}

extern "C"
void builtinHdgSourcewonly(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg,
                    const dstype* odg, const dstype* wdg, const dstype* uinf,
                    const dstype* param, dstype time, int builtinmodelID, int ng, int nc,
                    int ncu, int nd, int ncx, int nco, int ncw)
{
    switch (builtinmodelID) {
      case 1: m1::HdgSourcewonly(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 2: m2::HdgSourcewonly(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;     
      case 3: m3::HdgSourcewonly(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      case 4: m4::HdgSourcewonly(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, builtinmodelID, ng, nc, ncu, nd, ncx, nco, ncw); return;
      default: std::fprintf(stderr, "ERROR: Unknown builtinmodelID=%d in HdgSourcewonly\n", builtinmodelID); std::abort();
    }
}
