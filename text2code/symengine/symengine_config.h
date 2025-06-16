#ifndef SYMENGINE_CONFIG_HPP
#define SYMENGINE_CONFIG_HPP

#define SYMENGINE_MAJOR_VERSION 0
#define SYMENGINE_MINOR_VERSION 14
#define SYMENGINE_PATCH_VERSION 0
#define SYMENGINE_VERSION "0.14.0"

/* Define if you want to enable ASSERT testing in SymEngine */
/* #undef WITH_SYMENGINE_ASSERT */

/* Define if you want to enable GMP support in SymEngine */
/* #undef HAVE_SYMENGINE_GMP */

/* Define if you want to enable SYMENGINE_RCP support in SymEngine */
#define WITH_SYMENGINE_RCP

/* Define if you want to enable TEUCHOS support in SymEngine */
/* #undef WITH_SYMENGINE_TEUCHOS */

/* Define if you want to enable SYMENGINE_THREAD_SAFE support in SymEngine */
#define WITH_SYMENGINE_THREAD_SAFE

/* Define if you want to enable ECM support in SymEngine */
/* #undef HAVE_SYMENGINE_ECM */

/* Define if you want to enable PRIMESIEVE support in SymEngine */
/* #undef HAVE_SYMENGINE_PRIMESIEVE */

/* Define if you want to use virtual TypeIDs in SymEngine */
/* #undef WITH_SYMENGINE_VIRTUAL_TYPEID */

/* Define if you want to enable Flint support in SymEngine */
/* #undef HAVE_SYMENGINE_FLINT */

/* Define if you want to enable ARB support in SymEngine */
/* #undef HAVE_SYMENGINE_ARB */

/* Define if you want to enable MPFR support in SymEngine */
/* #undef HAVE_SYMENGINE_MPFR */

/* Define if you want to enable Piranha support in SymEngine */
/* #undef HAVE_SYMENGINE_PIRANHA */

/* Define if you want to enable BOOST support in SymEngine */
#define HAVE_SYMENGINE_BOOST

/* Define if you want to enable PTHREAD support in SymEngine */
/* #undef HAVE_SYMENGINE_PTHREAD */

/* Define if you want to enable MPC support in SymEngine */
/* #undef HAVE_SYMENGINE_MPC */

/* Define if you want to enable LLVM support in SymEngine */
/* #undef HAVE_SYMENGINE_LLVM */

/* Define if the C compiler supports __FUNCTION__ but not __func__ */
/* #undef HAVE_C_FUNCTION_NOT_FUNC */

/* Define if the C++ compiler supports default constructors */
#define HAVE_DEFAULT_CONSTRUCTORS

/* Define if the C++ compiler supports noexcept specifier */
#define HAVE_SYMENGINE_NOEXCEPT

/* Define if the C++ compiler supports std::is_constructible */
#define HAVE_SYMENGINE_IS_CONSTRUCTIBLE

/* Define if the C++ compiler supports std::unordered_map<>::reserve() */
#define HAVE_SYMENGINE_RESERVE

/* Define if the C++ compiler has std::to_string */
#define HAVE_SYMENGINE_STD_TO_STRING

/* Define if the C++ compiler has RTTI */
#define HAVE_SYMENGINE_RTTI 1

#define SYMENGINE_GMPXX 0
#define SYMENGINE_PIRANHA 1
#define SYMENGINE_FLINT 2
#define SYMENGINE_GMP 3
#define SYMENGINE_BOOSTMP 4

#define SYMENGINE_INTEGER_CLASS SYMENGINE_BOOSTMP

#define SYMENGINE_SIZEOF_LONG_DOUBLE 16

#ifdef HAVE_SYMENGINE_NOEXCEPT
#  define SYMENGINE_NOEXCEPT noexcept
#else
#  define SYMENGINE_NOEXCEPT
#endif

#include <symengine/symengine_export.h>

#ifdef __CLING__
#include "symengine/symengine_config_cling.h"
#endif

#endif
