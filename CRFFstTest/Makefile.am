AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFFstTest
CRFFstTest_SOURCES = src/Main.cpp
CRFFstTest_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
