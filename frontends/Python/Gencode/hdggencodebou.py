import sympy
from varsassign import varsassign
from sympyassign2 import sympyassign2

def hdggencodebou(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):

    cpufile = "Hdg" + filename
    tmp = "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    str_code = "\tKokkos::parallel_for(" + "\"" + filename + "\"" + ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += str_code

    fstr = str(f.flatten())
    str_code = ""

    # Call to the already implemented varsassign function
    str_code = varsassign(str_code, "param", len(param), 0, fstr)
    str_code = varsassign(str_code, "uinf", len(uinf), 0, fstr)
    str_code = varsassign(str_code, "tau", len(tau), 0, fstr)
    str_code = varsassign(str_code, "xdg", len(xdg), 1, fstr)
    str_code = varsassign(str_code, "udg", len(udg), 1, fstr)
    str_code = varsassign(str_code, "uhg", len(uhg), 1, fstr)
    str_code = varsassign(str_code, "odg", len(odg), 1, fstr)
    str_code = varsassign(str_code, "wdg", len(wdg), 1, fstr)
    str_code = varsassign(str_code, "nlg", len(nlg), 1, fstr)

    # Call to the already implemented sympyassign2 function
    str_code = sympyassign2(str_code, f.flatten(), udg, wdg, uhg)

    strkk += str_code + "\t});\n" + "}\n\n"
    strkk = strkk.replace("T ", "dstype ")

    return strkk
