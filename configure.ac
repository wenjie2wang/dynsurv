# dynsurv configure.ac by Jun Yan <jun.yan@uconn.edu>
# Adapted from 
# RQuantLib configure.in by Dirk Eddelbuettel <edd@debian.org>
#
# Using pieces borrowed from RPgSQL, GNU Gretl, GNU R and QuantLib
# Greatly simplified thanks to quantlib-config
# Another contribution by Kurt Hornik gratefully acknowledged

# require at least autoconf 2.50
AC_PREREQ(2.50)

# Process this file with autoconf to produce a configure script.
AC_INIT(dynsurv, 0.1-6)

# We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

AC_PROG_CXX
if test "${GXX}" = yes; then
    gxx_version=`${CXX} -v 2>&1 | grep "^.*g.. version" | \\
		       sed -e 's/^.*g.. version *//'`
    case ${gxx_version} in
        1.*|2.*)
	     AC_MSG_WARN([Only g++ version 3.0 or greater can be used with RQuantib.])
	     AC_MSG_ERROR([Please use a different compiler.])   
        ;;
	4.6.*|4.7.*|4.8.*)
	     gxx_newer_than_45="-fpermissive"
	;;
    esac
fi

AC_DEFUN(AC_PROG_R, [AC_CHECK_PROG(R,R,yes)])
AC_PROG_R

# Next segment by Kurt
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
CXX=`${R_HOME}/bin/R CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

## look for Rscript, but use the one found via R_HOME to allow for multiple installations
AC_DEFUN(AC_PROG_RSCRIPT, [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
AC_PROG_RSCRIPT

## Rcpp is not needed; it was for RQantLib
dnl #AC_MSG_CHECKING([for Rscript])
dnl if test x"${RSCRIPT}" == x"yes" ; then
dnl #     AC_MSG_RESULT([yes])    
dnl      ## Rcpp compile flag providing header directory containing Rcpp.h
dnl      #rcpp_cxxflags=`${R_HOME}/bin/Rscript -e 'Rcpp:::CxxFlags()'`
dnl      ## link flag providing libary as well as path to library, and optionally rpath
dnl      rcpp_ldflags=`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`
dnl      ## now use all these
dnl      #AC_SUBST([CXXFLAGS],["${CXXFLAGS} $rcpp_cxxflags"])
dnl      AC_SUBST([LDFLAGS],["${LDFLAGS} $rcpp_ldflags"])
dnl else
dnl #    AC_MSG_RESULT([no])    
dnl     echo "
dnl   Your installation does not appear to have Rscript installed.

dnl   Please make sure that you have a working and complete R installation.
dnl "
dnl     exit 1
dnl fi

## check for boost, this segment as well as segment below are 
## from QuantLib's configure.ac
## 
##
## Help the compiler find external software
AC_ARG_WITH([boost-include],
            AC_HELP_STRING([--with-boost-include=INCLUDE_PATH],
                           [Supply the location of Boost header files]),
            [boost_include_path=$withval],
            [boost_include_path=""])
if test [ -n "$boost_include_path" ] ; then
   AC_SUBST([CXXFLAGS],["${CXXFLAGS} -I${boost_include_path}"])
fi
AC_ARG_WITH([boost-lib],
            AC_HELP_STRING([--with-boost-lib=LIB_PATH],
                           [Supply the location of Boost libraries]),
            [boost_lib_path=$withval],
            [boost_lib_path=""])
if test [ -n "$boost_lib_path" ] ; then
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${boost_lib_path}"])
fi
##
## Check for Boost components
# CHECK_BOOST_DEVEL
# --------------------
# Check whether the Boost headers are available
AC_DEFUN([CHECK_BOOST_DEVEL],
[AC_MSG_CHECKING([for Boost development files])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>
     @%:@include <boost/shared_ptr.hpp>
     @%:@include <boost/assert.hpp>
     @%:@include <boost/current_function.hpp>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([Boost development files not found])
    ])
])

# CHECK_BOOST_VERSION
# ----------------------
# Check whether the Boost installation is up to date
AC_DEFUN([CHECK_BOOST_VERSION],
[AC_MSG_CHECKING([for minimal Boost version])
 AC_REQUIRE([CHECK_BOOST_DEVEL])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>],
    [@%:@if BOOST_VERSION < 103400
     @%:@error too old
     @%:@endif],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([outdated Boost installation])
    ])
])

# CHECK_BOOST
# ------------------------
# Boost-related tests
AC_DEFUN([CHECK_BOOST],
[AC_REQUIRE([CHECK_BOOST_DEVEL])
 AC_REQUIRE([CHECK_BOOST_VERSION])
])
CHECK_BOOST


AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_OUTPUT(src/Makevars)
echo "Completed configuration and ready to build."

### Local Variables: ***
### mode: autoconf ***
### autoconf-indentation: 4 ***
### End: ***


