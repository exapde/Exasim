#include "SymbolicFunctions.hpp"

std::vector<Expression> Flux(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> f;
    f.resize(3);

    Expression kappa = mu[0];
    f[0]  =  kappa*uq[1];
    f[1]  =  kappa*uq[2];
    f[2]  =  kappa*uq[3];
    return f;
}

std::vector<Expression> Source(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(1);

    s[0]  =  0.0;
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
    fb[0]  =  f[0]*n[0] + f[1]*n[1] + f[2]*n[2] + tau[0]*(uq[0]-uhat[0]);
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
    fb.resize(2);

    fb[0]  =  tau[0]*(mu[1] - uhat[0]);
    auto f = Flux(x, uq, v, w, eta, mu, t);
    fb[1]  =  f[0]*n[0] + f[1]*n[1] + f[2]*n[2] + tau[0]*(uq[0]-uhat[0]) - (mu[2]*n[0]+mu[3]*n[1]+mu[4]*n[2]);
    return fb;
}

std::vector<Expression> Initu(const std::vector<Expression>& x, const std::vector<Expression>& eta, const std::vector<Expression>& mu) {
    std::vector<Expression> ui;
    ui.resize(1);

    ui[0]  =  0.0;
    return ui;
}

