AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFLocalPosteriors
CRFLocalPosteriors_SOURCES = src/Main.cpp
CRFLocalPosteriors_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
