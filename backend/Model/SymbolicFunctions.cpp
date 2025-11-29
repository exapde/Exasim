#include "SymbolicFunctions.hpp"

std::vector<Expression> Flux(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> f;
    f.resize(2);

    Expression kappa = mu[0];
    f[0]  =  kappa*uq[1];
    f[1]  =  kappa*uq[2];
    return f;
}

std::vector<Expression> Source(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(1);

    Expression x1 = x[0];
    Expression x2 = x[1];
    s[0]  =  mul(mul(SymEngine::integer(2),mul(SymEngine::pi, SymEngine::pi)) , mul(sin(SymEngine::pi*x1),sin(SymEngine::pi*x2)));
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
    fb.resize(1);

    fb[0]  =  tau[0]*(0.0 - uhat[0]);
    return fb;
}

std::vector<Expression> Initu(const std::vector<Expression>& x, const std::vector<Expression>& eta, const std::vector<Expression>& mu) {
    std::vector<Expression> ui;
    ui.resize(1);

    ui[0]  =  0.0;
    return ui;
}

std::vector<Expression> VisScalars(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(2);

    s[0]  =  uq[0];
    s[1]  =  uq[1] + uq[2];
    return s;
}

std::vector<Expression> VisVectors(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(2);

    s[0]  =  uq[1];
    s[1]  =  uq[2];
    return s;
}

std::vector<Expression> QoIvolume(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> s;
    s.resize(2);

    Expression x1 = x[0];
    Expression x2 = x[1];
    auto t1 = Expression(SymEngine::pi);
    auto t2 = sin(t1*x1);
    auto t3 = sin(t1*x2);
    auto uexact = mul(t2,t3);
    s[0]  =  (uq[0] - uexact)*(uq[0] - uexact);
    s[1]  =  uq[0];
    return s;
}

std::vector<Expression> QoIboundary(const std::vector<Expression>& x, const std::vector<Expression>& uq, const std::vector<Expression>& v, const std::vector<Expression>& w, const std::vector<Expression>& uhat, const std::vector<Expression>& n, const std::vector<Expression>& tau, const std::vector<Expression>& eta, const std::vector<Expression>& mu, const Expression& t) {
    std::vector<Expression> fb;
    fb.resize(1);

    auto f = Flux(x, uq, v, w, eta, mu, t);
    fb[0]  =  f[0]*n[0] + f[1]*n[1] + tau[0]*(uq[0]-uhat[0]);
    return fb;
}

