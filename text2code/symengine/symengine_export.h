
#ifndef SYMENGINE_EXPORT_H
#define SYMENGINE_EXPORT_H

#ifdef SYMENGINE_STATIC_DEFINE
#  define SYMENGINE_EXPORT
#  define SYMENGINE_NO_EXPORT
#else
#  ifndef SYMENGINE_EXPORT
#    ifdef symengine_EXPORTS
        /* We are building this library */
#      define SYMENGINE_EXPORT 
#    else
        /* We are using this library */
#      define SYMENGINE_EXPORT 
#    endif
#  endif

#  ifndef SYMENGINE_NO_EXPORT
#    define SYMENGINE_NO_EXPORT 
#  endif
#endif

#ifndef SYMENGINE_DEPRECATED
#  define SYMENGINE_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SYMENGINE_DEPRECATED_EXPORT
#  define SYMENGINE_DEPRECATED_EXPORT SYMENGINE_EXPORT SYMENGINE_DEPRECATED
#endif

#ifndef SYMENGINE_DEPRECATED_NO_EXPORT
#  define SYMENGINE_DEPRECATED_NO_EXPORT SYMENGINE_NO_EXPORT SYMENGINE_DEPRECATED
#endif

/* NOLINTNEXTLINE(readability-avoid-unconditional-preprocessor-if) */
#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SYMENGINE_NO_DEPRECATED
#    define SYMENGINE_NO_DEPRECATED
#  endif
#endif

#endif /* SYMENGINE_EXPORT_H */
