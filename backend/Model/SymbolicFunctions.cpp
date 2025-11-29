#include "SymbolicFunctions.hpp"

std::vector<Expression> Flux(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> f;
    f.resize(8);

    Expression gam = mu[0];
    Expression gam1 = gam - 1.0;
    Expression Re = mu[1];
    Expression Pr = mu[2];
    Expression Minf = mu[3];
    Expression Re1 = 1/Re;
    Expression M2 = Minf * Minf;
    Expression c23 = 2.0/3.0;
    Expression fc = 1/(gam1*M2*Re*Pr);
    Expression r = uq[0];
    Expression ru = uq[1];
    Expression rv = uq[2];
    Expression rE = uq[3];
    Expression rx = uq[4];
    Expression rux = uq[5];
    Expression rvx = uq[6];
    Expression rEx = uq[7];
    Expression ry = uq[8];
    Expression ruy = uq[9];
    Expression rvy = uq[10];
    Expression rEy = uq[11];
    Expression r1 = 1/r;
    Expression uv = ru * r1;
    Expression vv = rv * r1;
    Expression E = rE * r1;
    Expression ke = 0.5*(uv*uv+vv*vv);
    Expression p = gam1*(rE-r*ke);
    Expression h = E+p*r1;
    Expression ux = (rux - rx*uv)*r1;
    Expression vx = (rvx - rx*vv)*r1;
    Expression kex = uv*ux + vv*vx;
    Expression px = gam1*(rEx - rx*ke - r*kex);
    Expression Tx = gam*M2*(px*r - p*rx)*r1*r1;
    Expression uy = (ruy - ry*uv)*r1;
    Expression vy = (rvy - ry*vv)*r1;
    Expression key = uv*uy + vv*vy;
    Expression py = gam1*(rEy - ry*ke - r*key);
    Expression Ty = gam*M2*(py*r - p*ry)*r1*r1;
    Expression txx = Re1*c23*(2*ux - vy);
    Expression txy = Re1*(uy + vx);
    Expression tyy = Re1*c23*(2*vy - ux);
    f[0]  =  ru;
    f[1]  =  ru*uv+p + txx;
    f[2]  =  rv*uv + txy;
    f[3]  =  ru*h + uv*txx + vv*txy + fc*Tx;
    f[4]  =  rv;
    f[5]  =  ru*vv + txy;
    f[6]  =  rv*vv+p + tyy;
    f[7]  =  rv*h + uv*txy + vv*tyy + fc*Ty;
    return f;
}

std::vector<Expression> Source(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(1);

    for (int i = 0; i < 4; ++i) {
         s[i] = Expression(0);
    }
    return s;
}

std::vector<Expression> Tdfunc(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> m;
    m.resize(1);

    for (int i = 0; i < 1; ++i) {
         m[i] = Expression(1);
    }
    return m;
}

std::vector<Expression> Fbou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> fb;
    fb.resize(1);

    auto f = Flux(x, uq, v, w, eta, mu, t);
    fb[0]  =  f[0]*n[0] + f[1]*n[1] + tau[0]*(uq[0]-uhat[0]);
    return fb;
}

std::vector<Expression> Ubou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> ub;
    ub.resize(1);

    ub[0]  =  0.0;
    return ub;
}

std::vector<Expression> FbouHdg(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> fb;
    fb.resize(8);

    auto f = Flux(x, uq, v, w, eta, mu, t);
    fb[0]  =  uq[0] - uhat[0];
    fb[1]  =  0.0  - uhat[1];
    fb[2]  =  0.0  - uhat[2];
    fb[3]  =  f[3]*n[0] + f[7]*n[1] + tau[0]*(uq[3]-uhat[3]);
    Expression gam = mu[0];
    Expression gam1 = gam - 1.0;
    Expression r = uhat[0];
    Expression ru = uhat[1];
    Expression rv = uhat[2];
    Expression rE = uhat[3];
    Expression nx = n[0];
    Expression ny = n[1];
    Expression r1 = 1/r;
    Expression uv = ru * r1;
    Expression vv = rv * r1;
    Expression E = rE * r1;
    Expression p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    Expression h = E+p*r1;
    auto a = sqrt(gam*p*r1);
    Expression run = ru*nx + rv*ny;
    Expression rut = -ru*ny + rv*nx;
    Expression un = run/r;
    Expression ut = rut/r;
    SymEngine::DenseMatrix K(4, 4);
    K.set(0, 0, Expression(1));
    K.set(1, 0, Expression(un-a));
    K.set(2, 0, Expression(ut));
    K.set(3, 0, Expression(h - un*a));
    K.set(0, 1, Expression(1));
    K.set(1, 1, Expression(un));
    K.set(2, 1, Expression(ut));
    K.set(3, 1, Expression(0.5*(un*un + ut*ut)));
    K.set(0, 2, Expression(0));
    K.set(1, 2, Expression(0));
    K.set(2, 2, Expression(1));
    K.set(3, 2, Expression(ut));
    K.set(0, 3, Expression(1));
    K.set(1, 3, Expression(un+a));
    K.set(2, 3, Expression(ut));
    K.set(3, 3, Expression(h+un*a));
    SymEngine::DenseMatrix Kinv(4, 4);
    K.inv(Kinv);
    SymEngine::DenseMatrix T(4, 4);
    T.set(0, 0, Expression(1));
    T.set(1, 0, Expression(0));
    T.set(2, 0, Expression(0));
    T.set(3, 0, Expression(0));
    T.set(0, 1, Expression(0));
    T.set(1, 1, Expression(nx));
    T.set(2, 1, Expression(-ny));
    T.set(3, 1, Expression(0));
    T.set(0, 2, Expression(0));
    T.set(1, 2, Expression(ny));
    T.set(2, 2, Expression(nx));
    T.set(3, 2, Expression(0));
    T.set(0, 3, Expression(0));
    T.set(1, 3, Expression(0));
    T.set(2, 3, Expression(0));
    T.set(3, 3, Expression(1));
    SymEngine::DenseMatrix Tinv(4, 4);
    T.inv(Tinv);
    SymEngine::DenseMatrix Lambda(4, 4);
    Lambda.set(0, 0, Expression(tanh(100*(un-a))));
    Lambda.set(1, 0, Expression(0));
    Lambda.set(2, 0, Expression(0));
    Lambda.set(3, 0, Expression(0));
    Lambda.set(0, 1, Expression(0));
    Lambda.set(1, 1, Expression(tanh(100*un)));
    Lambda.set(2, 1, Expression(0));
    Lambda.set(3, 1, Expression(0));
    Lambda.set(0, 2, Expression(0));
    Lambda.set(1, 2, Expression(0));
    Lambda.set(2, 2, Expression(tanh(100*un)));
    Lambda.set(3, 2, Expression(0));
    Lambda.set(0, 3, Expression(0));
    Lambda.set(1, 3, Expression(0));
    Lambda.set(2, 3, Expression(0));
    Lambda.set(3, 3, Expression(tanh(100*(un+a))));
    SymEngine::DenseMatrix L(4, 4);
    Tinv.mul_matrix(K, L);
    SymEngine::DenseMatrix R(4, 4);
    Kinv.mul_matrix(T, R);
    SymEngine::DenseMatrix Tmp(4, 4);
    L.mul_matrix(Lambda, Tmp);
    SymEngine::DenseMatrix An(4, 4);
    Tmp.mul_matrix(R, An);
    SymEngine::DenseMatrix bn(4, 1);
    for (int i = 0; i <= 3; ++i) {
    bn.set(i, 0, Expression(uq[i] - mu[4+i]));
    }
    SymEngine::DenseMatrix cn(4, 1);
    An.mul_matrix(bn, cn);
    std::vector<Expression> dn(4);
    for (int i = 0; i <= 3; ++i) {
    dn[i]  = cn.get(i, 0);
    }
    fb[4]  =  0.5*(uq[0] + mu[4] + dn[0]) - uhat[0];
    fb[5]  =  0.5*(uq[1] + mu[5] + dn[1]) - uhat[1];
    fb[6]  =  0.5*(uq[2] + mu[6] + dn[2]) - uhat[2];
    fb[7]  =  0.5*(uq[3] + mu[7] + dn[3]) - uhat[3];
    return fb;
}

std::vector<Expression> Initu(const std::vector<Expression>& x, const std::vector<Expression>& eta, const std::vector<Expression>& mu) {
    std::vector<Expression> ui;
    ui.resize(1);

    ui[0]  =  0.0;
    return ui;
}

