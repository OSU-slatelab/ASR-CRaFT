#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
LT_PREREQ([2.4.6])
AC_PREREQ([2.68])
AC_INIT([feacat], [1.1], [https://github.com/OSU-slatelab])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign])

#AC_CONFIG_SRCDIR([CRFTrain/src/Main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
gl_EARLY
AC_PROG_RANLIB
AM_PROG_AR
# libtool
# LT_INIT([dlopen])

# Checks for libraries.

# gnulib
gl_INIT
AC_CONFIG_MACRO_DIRS([m4])


AC_LIB_LINKFLAGS([quicknet3])
AC_LIB_LINKFLAGS([fst])

QN_LIBFILE=
QN_HEADERS=
AS_IF([test "x$ac_cv_libquicknet3_prefix" != "x"],
	    [QN_HEADERS="$LIBQUICKNET3_PREFIX/include/quicknet3";
	    QN_LIBFILE="$LIBQUICKNET3_PREFIX/lib/libquicknet3.a"],
	    AS_IF([test -f "/usr/lib/libquicknet3.a"],
 	    		[QN_HEADERS="/usr/include/quicknet3";
			QN_LIBFILE="/usr/lib/libquicknet3.a"],
			[QN_HEADERS="/usr/local/include/quicknet3";
			QN_LIBFILE="/usr/local/lib/libquicknet3.a"]))
AC_SUBST([QN_HEADERS])
AC_SUBST([QN_LIBFILE])

#echo "CPPFLAGS"
#echo $CPPFLAGS
#echo $INCQUICKNET3

# Checks for header files.
#AC_CHECK_HEADERS([malloc.h memory.h stdlib.h string.h sys/time.h])
#AC_CHECK_HEADERS([QuickNet.h])
# Checks for typedefs, structures, and compiler characteristics.
ACDPWE_CXX_BOOL
#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_SIZE_T
#AC_TYPE_UINT32_T

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
AC_CHECK_FUNCS([strdup])


# for libtool

#		 lib/Makefile

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
