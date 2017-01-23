bin_PROGRAMS = CRFTrain
CRFTrain_SOURCES = src/Main.cpp
CRFTrain_LDADD = $(top_builddir)/CRF/libCRF.a $(LIBQUICKNET3) $(LIBFST) -ldl -lpthread
CRFTrain_CPPFLAGS = -I$(top_srcdir)/CRF/src -I$(QN_HEADERS)
