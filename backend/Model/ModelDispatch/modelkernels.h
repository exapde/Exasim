/**
 * @file modelkernels.h
 * @brief Shared low-level kernel declarations and function-pointer types.
 *
 * This header declares the provider-side model kernel surface used by
 * ModelDriversDispatch and provider ABI tables. It intentionally stays at the
 * raw array/scalar argument level and does not depend on backend-facing driver
 * wrappers or solver runtime objects.
 */
#ifndef __EXASIM_MODEL_DISPATCH_MODELKERNELS_H__
#define __EXASIM_MODEL_DISPATCH_MODELKERNELS_H__

#include "../../Common/common.h"

void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ng,
                int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg,
                  const dstype* odg, const dstype* wdg, const dstype* uinf,
                  const dstype* param, dstype time, int modelnumber, int ng,
                  int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosSourcew(dstype* f, const dstype* xdg, const dstype* udg,
                   const dstype* odg, const dstype* wdg, const dstype* uinf,
                   const dstype* param, dstype time, int modelnumber, int ng,
                   int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                   int npe, int ne);

void KokkosTdfunc(dstype* f, const dstype* xdg, const dstype* udg,
                  const dstype* odg, const dstype* wdg, const dstype* uinf,
                  const dstype* param, dstype time, int modelnumber, int ng,
                  int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosAvfield(dstype* f, const dstype* xdg, const dstype* udg,
                   const dstype* odg, const dstype* wdg, const dstype* uinf,
                   const dstype* param, dstype time, int modelnumber, int ng,
                   int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                   int npe, int ne);

void KokkosEoS(dstype* f, const dstype* xdg, const dstype* udg,
               const dstype* odg, const dstype* wdg, const dstype* uinf,
               const dstype* param, dstype time, int modelnumber, int ng,
               int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
               int npe, int ne);

void KokkosEoSdu(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ng,
                 int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                 int npe, int ne);

void KokkosEoSdw(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ng,
                 int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                 int npe, int ne);

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uhg,
                const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ib,
                int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosUbou(dstype* f, const dstype* xdg, const dstype* udg,
                const dstype* odg, const dstype* wdg, const dstype* uhg,
                const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ib,
                int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosFhat(dstype* f, const dstype* xdg, const dstype* udg1,
                const dstype* udg2, const dstype* odg1, const dstype* odg2,
                const dstype* wdg1, const dstype* wdg2, const dstype* uhg,
                const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ng,
                int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosUhat(dstype* f, const dstype* xdg, const dstype* udg1,
                const dstype* udg2, const dstype* odg1, const dstype* odg2,
                const dstype* wdg1, const dstype* wdg2, const dstype* uhg,
                const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ng,
                int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosStab(dstype* f, const dstype* xdg, const dstype* udg1,
                const dstype* udg2, const dstype* odg1, const dstype* odg2,
                const dstype* wdg1, const dstype* wdg2, const dstype* uhg,
                const dstype* nlg, const dstype* tau, const dstype* uinf,
                const dstype* param, dstype time, int modelnumber, int ng,
                int nc, int ncu, int nd, int ncx, int nco, int ncw);

void KokkosOutput(dstype* f, const dstype* xdg, const dstype* udg,
                  const dstype* odg, const dstype* wdg, const dstype* uinf,
                  const dstype* param, dstype time, int modelnumber, int ng,
                  int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                  int npe, int ne);

void KokkosMonitor(dstype* f, const dstype* xdg, const dstype* udg,
                   const dstype* odg, const dstype* wdg, const dstype* uinf,
                   const dstype* param, dstype time, int modelnumber, int ng,
                   int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                   int npe, int ne);

void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param, dstype time,
                      int modelnumber, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw);

void KokkosVisVectors(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param, dstype time,
                      int modelnumber, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw);

void KokkosVisTensors(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param, dstype time,
                      int modelnumber, int ng, int nc, int ncu, int nd,
                      int ncx, int nco, int ncw);

void KokkosQoIvolume(dstype* f, const dstype* xdg, const dstype* udg,
                     const dstype* odg, const dstype* wdg,
                     const dstype* uinf, const dstype* param, dstype time,
                     int modelnumber, int ng, int nc, int ncu, int nd,
                     int ncx, int nco, int ncw);

void KokkosQoIboundary(dstype* f, const dstype* xdg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nlg,
                       const dstype* tau, const dstype* uinf,
                       const dstype* param, dstype time, int modelnumber,
                       int ib, int ng, int nc, int ncu, int nd, int ncx,
                       int nco, int ncw);

void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf,
                 const dstype* param, int modelnumber, int ng, int ncx,
                 int nfield, int npe, int ne);

void KokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf,
                 const dstype* param, int modelnumber, int ng, int ncx,
                 int nfield, int npe, int ne);

void KokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf,
                   const dstype* param, int modelnumber, int ng, int ncx,
                   int nfield, int npe, int ne);

void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf,
                   const dstype* param, int modelnumber, int ng, int ncx,
                   int nfield, int npe, int ne);

void KokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf,
                   const dstype* param, int modelnumber, int ng, int ncx,
                   int nfield, int npe, int ne);

void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf,
              const dstype* param, int modelnumber, int ng, int ncx,
              int nfield_runtime, int npe, int ne);

void cpuInitq(dstype* f, const dstype* xdg, const dstype* uinf,
              const dstype* param, int modelnumber, int ng, int ncx,
              int nfield_runtime, int npe, int ne);

void cpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf,
                const dstype* param, int modelnumber, int ng, int ncx,
                int nfield_runtime, int npe, int ne);

void cpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf,
                const dstype* param, int modelnumber, int ng, int ncx,
                int nfield_runtime, int npe, int ne);

void cpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf,
                const dstype* param, int modelnumber, int ng, int ncx,
                int nfield_runtime, int npe, int ne);

void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg,
             const dstype* udg, const dstype* odg, const dstype* wdg,
             const dstype* uinf, const dstype* param, dstype time,
             int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
             int nco, int ncw);

void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg,
               const dstype* udg, const dstype* odg, const dstype* wdg,
               const dstype* uinf, const dstype* param, dstype time,
               int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
               int nco, int ncw);

void HdgSourcew(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg,
                const dstype* udg, const dstype* odg, const dstype* wdg,
                const dstype* uinf, const dstype* param, dstype time,
                int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
                int nco, int ncw);

void HdgSourcewonly(dstype* f, dstype* f_wdg, const dstype* xdg,
                    const dstype* udg, const dstype* odg, const dstype* wdg,
                    const dstype* uinf, const dstype* param, dstype time,
                    int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
                    int nco, int ncw);

void HdgEoS(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg,
            const dstype* udg, const dstype* odg, const dstype* wdg,
            const dstype* uinf, const dstype* param, dstype time,
            int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
            int nco, int ncw);

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
             const dstype* xdg, const dstype* udg, const dstype* odg,
             const dstype* wdg, const dstype* uhg, const dstype* nlg,
             const dstype* tau, const dstype* uinf, const dstype* param,
             dstype time, int modelnumber, int ib, int ng, int nc, int ncu,
             int nd, int ncx, int nco, int ncw);

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ib,
                 int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

void HdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
             const dstype* xdg, const dstype* udg, const dstype* odg,
             const dstype* wdg, const dstype* uhg, const dstype* nlg,
             const dstype* tau, const dstype* uinf, const dstype* param,
             dstype time, int modelnumber, int ib, int ng, int nc, int ncu,
             int nd, int ncx, int nco, int ncw);

void HdgFintonly(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ib,
                 int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

void HdgFext(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
             const dstype* xdg, const dstype* udg, const dstype* odg,
             const dstype* wdg, const dstype* uhg, const dstype* nlg,
             const dstype* uext, const dstype* tau, const dstype* uinf,
             const dstype* param, dstype time, int modelnumber, int ib, int ng,
             int nc, int ncu, int nd, int ncx, int nco, int ncw);

void HdgFextonly(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* uext, const dstype* tau,
                 const dstype* uinf, const dstype* param, dstype time,
                 int modelnumber, int ib, int ng, int nc, int ncu, int nd,
                 int ncx, int nco, int ncw);

#endif
