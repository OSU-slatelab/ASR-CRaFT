bin_PROGRAMS = CRFFstDecode
CRFFstDecode_SOURCES = src/Main.cpp
CRFFstDecode_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl
CRFFstDecode_CPPFLAGS = -I$(top_srcdir)/CRF/src -I$(QN_HEADERS)
