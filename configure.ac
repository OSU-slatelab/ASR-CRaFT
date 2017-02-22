#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
LT_PREREQ([2.4.6])
AC_PREREQ([2.68])
AC_INIT([ASR-CRaFT], [1.1], [https://github.com/OSU-slatelab])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign])

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
# check $ac_cv_libquicknet3_prefix; if not set, libquicknet3
# was found in /usr/lib or /usr/local/lib (latter of which is its
# default install location). Set QN_HEADERS accordingly. Don't want
# to source QNConfig.sh since it defines $prefix, which might hose
# the configure script. This whole thing seems hacky, but I'm just
# working with what quicknet gives me.

# also, we need the actual file of the quicknet library for linking,
# since ar doesn't seem to like if automake passes it linker flags.

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

#echo $LIBQUICKNET3
#echo $LIBQUICKNET3_PREFIX

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
