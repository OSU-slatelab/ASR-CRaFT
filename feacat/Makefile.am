bin_PROGRAMS = feacat labcat
feacat_SOURCES = feacat.C Range.C cle.c cle.h Range.H
feacat_LDADD = $(LIBQUICKNET3)
feacat_CPPFLAGS = -I$(QN_HEADERS)
labcat_SOURCES = feacat.C Range.C cle.c
labcat_LDADD = $(LIBQUICKNET3)
labcat_CPPFLAGS = -I$(QN_HEADERS) -DLABCAT
dist_man1_MANS = feacat.man labcat.man
