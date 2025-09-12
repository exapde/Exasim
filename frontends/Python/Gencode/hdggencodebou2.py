import sympy
from varsassign import varsassign
from getccode import getccode

def hdggencodebou2(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):
    cpufile = "Hdg" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
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

    # Call to the already implemented getccode function
    str1 = getccode(f.flatten(), "f")
    str_code += str1

    strkk += str_code + "\t});\n" + "}\n\n"
    strkk = strkk.replace("T ", "dstype ")

    return strkk

# Example usage:
# f = sympy.symbols('f0:10')
# xdg = sympy.symbols('xdg0:10')
# udg = sympy.symbols('udg0:10')
# odg = sympy.symbols('odg0:10')
# wdg = sympy.symbols('wdg0:10')
# uhg = sympy.symbols('uhg0:10')
# nlg = sympy.symbols('nlg0:10')
# tau = sympy.symbols('tau0:10')
# uinf = sympy.symbols('uinf0:10')
# param = sympy.symbols('param0:10')
# time = sympy.symbols('time')
# filename = "example_filename"
# code = hdggencodebou2(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)
# print(code)
