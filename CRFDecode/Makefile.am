AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFDecode
CRFDecode_SOURCES = src/Main.cpp 
CRFDecode_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
