bin_PROGRAMS = CRFDecode
CRFDecode_SOURCES = src/Main.cpp 
CRFDecode_LDADD = $(top_builddir)/CRF/libCRF.a $(LIBQUICKNET3) $(LIBFST) -ldl
CRFDecode_CPPFLAGS = -I$(top_srcdir)/CRF/src -I$(QN_HEADERS)
