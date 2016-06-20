AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFGatherStats
CRFGatherStats_SOURCES = src/Main.cpp
CRFGatherStats_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
