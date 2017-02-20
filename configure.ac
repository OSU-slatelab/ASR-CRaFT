#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
LT_PREREQ([2.4.6])
AC_PREREQ([2.68])
AC_INIT([ASR-CRaFT], [1.1], [https://github.com/OSU-slatelab])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign])

#AC_CONFIG_SRCDIR([CRFTrain/src/Main.cpp])
#AC_CONFIG_HEADERS([config.h])

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
#AC_CONFIG_MACRO_DIRS([m4])

AC_LIB_LINKFLAGS([quicknet3])
AC_LIB_LINKFLAGS([fst])
QN_HEADERS=`echo $LIBQUICKNET3 | sed 's_lib/libquicknet3.a$_include/quicknet3_g'`
AC_SUBST([QN_HEADERS])

#echo "CPPFLAGS"
#echo $CPPFLAGS
#echo $INCQUICKNET3

AC_ARG_WITH([kaldi-timit-s5],
	    AS_HELP_STRING([--with-kaldi-timit-s5=DIR],
	    [optional directory location of an installed kaldi timit s5 experiment, for running the demo]),
	    [],
	    [with_kaldi_timit_s5=no])
TIMITDIR=
AS_IF([test "x$with_kaldi_timit_s5" != xno],
	    [AC_SUBST([TIMITDIR], [$with_kaldi_timit_s5])])

# Checks for header files.
AC_CHECK_HEADERS([malloc.h memory.h stdlib.h string.h])
#AC_CHECK_HEADERS([QuickNet.h])
# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_SIZE_T
#AC_TYPE_UINT32_T

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([gettimeofday memset sqrt])


# for libtool

#		 lib/Makefile
AC_CONFIG_SUBDIRS([feacat])
AC_CONFIG_FILES([CRF/Makefile
		 CRFDecode/Makefile
                 CRFFstDecode/Makefile
                 CRFTrain/Makefile
		 demo/Makefile
		 demo/kaldi-mods/Makefile
		 Makefile])
AC_OUTPUT
