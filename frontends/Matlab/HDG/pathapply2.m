function uh = pathapply2(A, B1, B2, C1, C2, D1, D2, DL, DU, b, fpath, fintf, nep)

uh = 0*b;
[uh1, uh2] = pathextract2(b, fpath, fintf, nep);
[uh1, uh2] = pathsolve2(A, B1, B2, C1, C2, D1, D2, DL, DU, uh1, uh2);
uh = pathinsert2(uh, uh1, uh2, fpath, fintf, nep);

