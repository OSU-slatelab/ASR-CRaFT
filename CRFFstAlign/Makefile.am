AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFFstAlign
CRFFstAlign_SOURCES = src/Main.cpp
CRFFstAlign_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
