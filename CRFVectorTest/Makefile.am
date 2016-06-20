AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFVectorTest
CRFVectorTest_SOURCES = src/Main.cpp
CRFVectorTest_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
