AM_CPPFLAGS = "-I$(top_srcdir)/CRF/src"
bin_PROGRAMS = CRFTrain
CRFTrain_SOURCES = src/Main.cpp
CRFTrain_LDADD = $(top_builddir)/CRF/libCRF.a -lquicknet3 -lfst -ldl 
