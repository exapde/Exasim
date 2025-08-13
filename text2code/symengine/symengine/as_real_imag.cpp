#include <symengine/visitor.h>
#include <symengine/basic.h>

namespace SymEngine
{

void pow_number(const RCP<const Basic> &in_re, const RCP<const Basic> &in_im,
                unsigned long n, Ptr<RCP<const Basic>> &out_re,
                Ptr<RCP<const Basic>> &out_im)
{
    unsigned long mask = 1;
    RCP<const Basic> tmp;
    RCP<const Basic> p_re = in_re;
    RCP<const Basic> p_im = in_im;
    *out_re = one;
    *out_im = zero;

    while (true) {
        if (n & mask) {
            // Multiply r by p
            tmp = sub(mul(*out_re, p_re), mul(*out_im, p_im));
            *out_im = add(mul(*out_re, p_im), mul(*out_im, p_re));
            *out_re = tmp;
        }
        mask = mask << 1;
        if (mask > 0 and n >= mask) {
            // Multiply p by p
            tmp = sub(mul(p_re, p_re), mul(p_im, p_im));
            p_im = mul(two, mul(p_re, p_im));
            p_re = tmp;
        } else {
            break;
        }
    }
}

class RealImagVisitor : public BaseVisitor<RealImagVisitor>
{
private:
    Ptr<RCP<const Basic>> real_, imag_;

public:
    RealImagVisitor(const Ptr<RCP<const Basic>> &real,
                    const Ptr<RCP<const Basic>> &imag)
        : real_{real}, imag_{imag}
    {
    }

    void apply(const Basic &b)
    {
        b.accept(*this);
    }

    void bvisit(const Mul &x)
    {
        RCP<const Basic> rest = one;
        RCP<const Basic> fre_ = one, fim_ = zero;

        for (const auto &arg : x.get_args()) {
            apply(*arg);
            std::tie(fre_, fim_)
                = std::make_tuple(sub(mul(fre_, *real_), mul(fim_, *imag_)),
                                  add(mul(fre_, *imag_), mul(fim_, *real_)));
        }
        *real_ = fre_;
        *imag_ = fim_;
    }

    void bvisit(const Add &x)
    {
        RCP<const Basic> t;
        umap_basic_num dr, dim;
        RCP<const Number> coefr = zero, coefim = zero, coef;

        for (const auto &arg : x.get_args()) {
            apply(*arg);
            if (is_a_Number(**real_)) {
                iaddnum(outArg(coefr), rcp_static_cast<const Number>(*real_));
            } else {
                Add::as_coef_term(*real_, outArg(coef), outArg(t));
                Add::dict_add_term(dr, coef, t);
            }
            if (is_a_Number(**imag_)) {
                iaddnum(outArg(coefim), rcp_static_cast<const Number>(*imag_));
            } else {
                Add::as_coef_term(*imag_, outArg(coef), outArg(t));
                Add::dict_add_term(dim, coef, t);
            }
        }

        *real_ = Add::from_dict(coefr, std::move(dr));
        *imag_ = Add::from_dict(coefim, std::move(dim));
    }

    void bvisit(const Pow &x)
    {
        RCP<const Basic> exp_;
        exp_ = x.get_exp();
        apply(*x.get_base());

        if (eq(**imag_, *zero)) {
            *real_ = x.rcp_from_this();
            *imag_ = zero;
            return;
        }
        if (is_a<Integer>(*exp_)) {
            if (static_cast<const Integer &>(*exp_).is_negative()) {
                auto magn = add(mul(*real_, *real_), mul(*imag_, *imag_));
                *imag_ = neg(*imag_);
                RCP<const Integer> expx = rcp_static_cast<const Integer>(exp_);
                expx = static_cast<const Integer &>(*exp_).neg();
                unsigned long n = numeric_cast<unsigned long>(
                    mp_get_ui(expx->as_integer_class()));
                RCP<const Basic> real1 = *real_, imag1 = *imag_;
                pow_number(real1, imag1, n, real_, imag_);
                magn = pow(magn, expx);
                *real_ = div(*real_, magn);
                *imag_ = div(*imag_, magn);
            } else {
                RCP<const Integer> expx = rcp_static_cast<const Integer>(exp_);
                unsigned long n = numeric_cast<unsigned long>(
                    mp_get_ui(expx->as_integer_class()));
                RCP<const Basic> real1 = *real_, imag1 = *imag_;
                pow_number(real1, imag1, n, real_, imag_);
            }
        } else if (is_a<Rational>(*exp_)) {
            auto magn = sqrt(add(mul(*real_, *real_), mul(*imag_, *imag_)));
            auto ang = atan2(*imag_, *real_);
            magn = pow(magn, exp_);
            ang = mul(ang, exp_);
            *real_ = mul(magn, cos(ang));
            *imag_ = mul(magn, sin(ang));
        } else {
            throw SymEngineException("Not Implemented");
        }
    }

    void bvisit(const ComplexBase &x)
    {
        *real_ = x.real_part();
        *imag_ = x.imaginary_part();
    }

    void bvisit(const Infty &x)
    {
        if (eq(x, *ComplexInf)) {
            *real_ = Nan;
            *imag_ = Nan;
        } else {
            *real_ = x.rcp_from_this();
            *imag_ = zero;
        }
    }

    void bvisit(const Sin &x)
    {
        apply(*x.get_arg());
        std::tie(*real_, *imag_) = std::make_tuple(
            mul(sin(*real_), cosh(*imag_)), mul(sinh(*imag_), cos(*real_)));
    }

    void bvisit(const Cos &x)
    {
        apply(*x.get_arg());
        std::tie(*real_, *imag_)
            = std::make_tuple(mul(cos(*real_), cosh(*imag_)),
                              neg(mul(sinh(*imag_), sin(*real_))));
    }

    void bvisit(const Tan &x)
    {
        apply(*x.get_arg());
        if (eq(**imag_, *zero)) {
            *real_ = x.rcp_from_this();
            *imag_ = zero;
            return;
        }
        auto twice_real_ = mul(two, *real_), twice_imag_ = mul(two, *imag_);
        auto den = add(cos(twice_real_), cosh(twice_imag_));
        *real_ = div(sin(twice_real_), den);
        *imag_ = div(sinh(twice_imag_), den);
    }

    void bvisit(const Csc &x)
    {
        apply(*div(one, sin(x.get_arg())));
    }

    void bvisit(const Sec &x)
    {
        apply(*div(one, cos(x.get_arg())));
    }

    void bvisit(const Cot &x)
    {
        apply(*x.get_arg());
        if (eq(**imag_, *zero)) {
            *real_ = x.rcp_from_this();
            return;
        }
        auto twice_real_ = mul(two, *real_), twice_imag_ = mul(two, *imag_);
        auto den = sub(cos(twice_real_), cosh(twice_imag_));
        *real_ = neg(div(sin(twice_real_), den));
        *imag_ = neg(div(sinh(twice_imag_), den));
    }

    void bvisit(const Sinh &x)
    {
        apply(*x.get_arg());
        std::tie(*real_, *imag_) = std::make_tuple(
            mul(sinh(*real_), cos(*imag_)), mul(sin(*imag_), cosh(*real_)));
    }

    void bvisit(const Cosh &x)
    {
        apply(*x.get_arg());
        std::tie(*real_, *imag_) = std::make_tuple(
            mul(cosh(*real_), cos(*imag_)), mul(sin(*imag_), sinh(*real_)));
    }

    void bvisit(const Tanh &x)
    {
        apply(*x.get_arg());
        if (eq(**imag_, *zero)) {
            *real_ = x.rcp_from_this();
            return;
        }
        auto sinh_re = sinh(*real_), cos_im = cos(*imag_);
        auto den = add(pow(sinh_re, two), pow(cos_im, two));
        *real_ = div(mul(sinh_re, cosh(*real_)), den);
        *imag_ = div(mul(sin(*imag_), cos_im), den);
    }

    void bvisit(const Csch &x)
    {
        apply(*div(one, sinh(x.get_arg())));
    }

    void bvisit(const Sech &x)
    {
        apply(*div(one, cosh(x.get_arg())));
    }

    void bvisit(const Coth &x)
    {
        apply(*x.get_arg());
        if (eq(**imag_, *zero)) {
            *real_ = x.rcp_from_this();
            return;
        }
        auto sinh_re = sinh(*real_), sin_im = sin(*imag_);
        auto den = add(pow(sinh_re, two), pow(sin_im, two));
        *real_ = div(mul(sinh_re, cosh(*real_)), den);
        *imag_ = neg(div(mul(sin_im, cos(*imag_)), den));
    }

    void bvisit(const Abs &x)
    {
        *real_ = x.rcp_from_this();
        *imag_ = zero;
    }

    void bvisit(const Function &x)
    {
        throw SymEngineException(
            "Not Implemented classes for real and imaginary parts");
    }

    void bvisit(const Symbol &x)
    {
        throw SymEngineException(
            "Not Implemented classes for real and imaginary parts");
    }

    void bvisit(const Basic &x)
    {
        *real_ = x.rcp_from_this();
        *imag_ = zero;
    }
};

void as_real_imag(const RCP<const Basic> &x, const Ptr<RCP<const Basic>> &real,
                  const Ptr<RCP<const Basic>> &imag)
{
    RealImagVisitor v(real, imag);
    v.apply(*x);
}

} // namespace SymEngine
