#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ASR-CRaFT], [1.1], [https://github.com/OSU-slatelab])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#AC_CONFIG_SRCDIR([CRFTrain/src/Main.cpp])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CONFIG_FILES([
	Makefile
	CRF/Makefile
	CRFDecode/Makefile
	CRFFstAlign/Makefile
	CRFFstDecode/Makefile
	CRFFstTest/Makefile
	CRFGatherStats/Makefile
	CRFLocalPosteriors/Makefile
	CRFLogMathTest/Makefile
	CRFTrain/Makefile
	CRFVectorTest/Makefile
])

# Checks for libraries.

# Checks for header files.
# AC_CHECK_HEADERS([malloc.h memory.h stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
# AC_CHECK_HEADER_STDBOOL
# AC_C_INLINE
# AC_TYPE_SIZE_T
# AC_TYPE_UINT32_T

# Checks for library functions.
# AC_FUNC_ERROR_AT_LINE
# AC_FUNC_MALLOC
# AC_CHECK_FUNCS([gettimeofday memset sqrt])

AC_OUTPUT
