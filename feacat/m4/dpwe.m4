dnl dpwe''s aclocal.m4     -*- sh -*-
dnl
dnl A place for my own autoconf macros
dnl 
dnl 1998jan23 Dan Ellis dpwe@icsi.berkeley.edu
dnl $Header: /u/drspeech/repos/feacat/aclocal.m4,v 1.2 2004/07/24 02:17:37 davidj Exp $
dnl
dnl Clients:
dnl  dpwetcl dpweutils_tcl farray_otcl gdtcl libdat otcl pfif_otcl
dnl  sound_otcl tclsh-readline aprl feacat
#--------------------------------------------------------------------
# ACDPWE_CXX_BOOL:
#       If the C++ compiler knows the bool type, define HAVE_BOOL.
#       Else don't.  This one uses the previously-defined $CXX
#--------------------------------------------------------------------

AC_DEFUN([ACDPWE_CXX_BOOL], [
AC_MSG_CHECKING([if C++ compiler supports "bool"])
cat > conftest.cc <<EOF
[#]line __oline__ "configure"
#include "confdefs.h"
int main() {
bool b=0;
return 0;
}
EOF
ac_try="$CXX -o conftest conftest.cc >/dev/null 2>conftest.out"
AC_TRY_EVAL([ac_try])
ac_err=`grep -v '^ *+' conftest.out`
if test -z "$ac_err"; then
  AC_DEFINE([HAVE_BOOL],[1], [Define if C++ compiler supports bool])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no]);
  echo "$ac_err" >&AC_FD_CC
fi
rm -f conftest*
	 ])




