#include "SymbolicFunctions.hpp"

std::vector<Expression> Flux(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> f;
    f.resize(8);

    Expression gam = mu[0];
    Expression gam1 = gam - 1.0;
    Expression Re = mu[1];
    Expression Pr = mu[2];
    Expression Minf = mu[3];
    Expression Tref = mu[9];
    Expression muRef = 1/Re;
    Expression M2 = Minf * Minf;
    Expression c23 = 2.0/3.0;
    Expression pinf = 1.0/(gam*M2);
    Expression Tinf = pinf/(gam-1.0);
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
    Expression av = v[0];
    Expression r1 = 1/r;
    Expression uv = ru * r1;
    Expression vv = rv * r1;
    Expression E = rE * r1;
    Expression ke = 0.5*(uv*uv+vv*vv);
    Expression p = gam1*(rE-r*ke);
    Expression h = E+p*r1;
    Expression T = p/(gam1*r);
    Expression Tphys = Tref/Tinf * T;
    Expression Ts = 110.4;
    Expression Tr = Tphys/Tref;
    Expression muphys = muRef * sqrt(Tr*Tr*Tr) * (Tref + Ts)/(Tphys + Ts);
    Expression fc = muphys*gam/Pr;
    Expression ux = (rux - rx*uv)*r1;
    Expression vx = (rvx - rx*vv)*r1;
    Expression kex = uv*ux + vv*vx;
    Expression px = gam1*(rEx - rx*ke - r*kex);
    Expression Tx = (px*r - p*rx)*r1*r1/gam1;
    Expression uy = (ruy - ry*uv)*r1;
    Expression vy = (rvy - ry*vv)*r1;
    Expression key = uv*uy + vv*vy;
    Expression py = gam1*(rEy - ry*ke - r*key);
    Expression Ty = (py*r - p*ry)*r1*r1/gam1;
    Expression txx = muphys*c23*(2*ux - vy);
    Expression txy = muphys*(uy + vx);
    Expression tyy = muphys*c23*(2*vy - ux);
    f[0]  =  ru + av*rx;
    f[1]  =  ru*uv+p + txx + av*rux;
    f[2]  =  rv*uv + txy + av*rvx;
    f[3]  =  ru*h + uv*txx + vv*txy + fc*Tx + av*rEx;
    f[4]  =  rv + av*ry;
    f[5]  =  ru*vv + txy + av*ruy;
    f[6]  =  rv*vv+p + tyy + av*rvy;
    f[7]  =  rv*h + uv*txy + vv*tyy + fc*Ty + av*rEy;
    return f;
}

std::vector<Expression> Source(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(4);

    s[0]  =  0;
    s[1]  =  0;
    s[2]  =  0;
    s[3]  =  0;
    return s;
}

std::vector<Expression> Tdfunc(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> m;
    m.resize(4);

    for (int i = 0; i < 4; ++i) {
         m[i] = Expression(1);
    }
    return m;
}

std::vector<Expression> Fbou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> fb;
    fb.resize(8);

    for (int i = 0; i < 8; ++i) {
         fb[i] = Expression(0);
    }
    return fb;
}

std::vector<Expression> Ubou(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> ub;
    ub.resize(8);

    for (int i = 0; i < 8; ++i) {
         ub[i] = Expression(0);
    }
    return ub;
}

std::vector<Expression> FbouHdg(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> fb;
    fb.resize(12);

    Expression gam = mu[0];
    Expression gam1 = gam - 1.0;
    Expression Tinf = mu[8];
    Expression Tref = mu[9];
    Expression Twall = mu[10];
    Expression TisoW = Twall/Tref * Tinf;
    fb[0]  =  mu[4] - uhat[0];
    fb[1]  =  mu[5] - uhat[1];
    fb[2]  =  mu[6] - uhat[2];
    fb[3]  =  mu[7] - uhat[3];
    fb[4]  =  uq[0] - uhat[0];
    fb[5]  =  uq[1] - uhat[1];
    fb[6]  =  uq[2] - uhat[2];
    fb[7]  =  uq[3] - uhat[3];
    fb[8]  =   uq[0] - uhat[0];
    fb[9]  =   0.0  - uhat[1];
    fb[10]  =  0.0  - uhat[2];
    fb[11]  =  uhat[0]*TisoW - uhat[3];
    return fb;
}

std::vector<Expression> Initu(const std::vector<Expression>& x, const std::vector<Expression>& eta, const std::vector<Expression>& mu) {
    std::vector<Expression> ui;
    ui.resize(4);

    ui[0]  =  mu[4];
    ui[1]  =  mu[5];
    ui[2]  =  mu[6];
    ui[3]  =  mu[7];
    return ui;
}

std::vector<Expression> VisScalars(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(4);

    s[0]  =  uq[0];
    s[1]  =  uq[1]/uq[0];
    s[2]  =  uq[2]/uq[0];
    s[3]  =  0.4*(uq[3] - 0.5*(uq[1]*s[1] + uq[2]*s[2]));
    return s;
}

std::vector<Expression> VisVectors(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(2);

    s[0]  =  uq[1];
    s[1]  =  uq[2];
    return s;
}

