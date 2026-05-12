void HdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
    (void)f; (void)f_udg; (void)f_wdg; (void)f_uhg; (void)xdg; (void)udg;
    (void)odg; (void)wdg; (void)uhg; (void)nlg; (void)tau; (void)uinf;
    (void)param; (void)time; (void)modelnumber; (void)ib; (void)ng; (void)nc;
    (void)ncu; (void)nd; (void)ncx; (void)nco; (void)ncw;
    printf("UserDefined HdgFint is not implemented. Multi-domain coupling remains deferred.\n");
    exit(-1);
}
