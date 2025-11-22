#pragma once
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/integer.h>

using SymEngine::Expression;
using SymEngine::integer;


// Wrapper functions to allow the user to use math functions in the global namespace while mapping them to SymEngine functions that return Expression objects.
// Many SymEngine overloads are only defined for SymEngine::Expression objects, but the SymEngine math functions return RCP objects, so they are incompatible with other SymEngine::Expressions in following expressions. These wrappers solve this problem.

static inline Expression sin(const Expression &c) {
    return Expression(SymEngine::sin(c));
}

static inline Expression cos(const Expression &c) {
    return Expression(SymEngine::cos(c));
}

static inline Expression tan(const Expression &c) {
    return Expression(SymEngine::tan(c));
}

static inline Expression asin(const Expression &c) {
    return Expression(SymEngine::asin(c));
}

static inline Expression acos(const Expression &c) {
    return Expression(SymEngine::acos(c));
}

static inline Expression atan(const Expression &c) {
    return Expression(SymEngine::atan(c));
}

static inline Expression atan2(const Expression &y, const Expression &x) {
    return Expression(SymEngine::atan2(y, x));
}

static inline Expression sinh(const Expression &c) {
    return Expression(SymEngine::sinh(c));
}

static inline Expression cosh(const Expression &c) {
    return Expression(SymEngine::cosh(c));
}

static inline Expression tanh(const Expression &c) {
    return Expression(SymEngine::tanh(c));
}

static inline Expression asinh(const Expression &c) {
    return Expression(SymEngine::asinh(c));
}

static inline Expression acosh(const Expression &c) {
    return Expression(SymEngine::acosh(c));
}

static inline Expression atanh(const Expression &c) {
    return Expression(SymEngine::atanh(c));
}

static inline Expression exp(const Expression &c) {
    return Expression(SymEngine::exp(c));
}

static inline Expression log(const Expression &c) {
    return Expression(SymEngine::log(c));
}

static inline Expression log10(const Expression &c) {
    return Expression(SymEngine::log(c)) / Expression(SymEngine::log(integer(10)));
}

static inline Expression pow(const Expression &base, const Expression &exp) {
    return Expression(SymEngine::pow(base, exp));
}

static inline Expression sqrt(const Expression &c) {
    return Expression(SymEngine::sqrt(c));
}

static inline Expression abs(const Expression &c) {
    return Expression(SymEngine::abs(c));
}

static inline Expression floor(const Expression &c) {
    return Expression(SymEngine::floor(c));
}

static inline Expression ceiling(const Expression &c) {
    return Expression(SymEngine::ceiling(c));
}

static inline Expression erf(const Expression &c) {
    return Expression(SymEngine::erf(c));
}

static inline Expression erfc(const Expression &c) {
    return Expression(SymEngine::erfc(c));
}

static inline Expression gamma(const Expression &c) {
    return Expression(SymEngine::gamma(c));
}

static inline Expression lgamma(const Expression &c) {
    return Expression(SymEngine::loggamma(c));
}
