##
## Interface: ngspice 1.0.0.0
## created by: admsXml-2.3.0 - Sun, 21 Sep 2014 18:23:46

## Process this file with automake to produce Makefile.in

ADMSXMLINTERFACE=$(srcdir)/../admst

noinst_LTLIBRARIES = libhicum2.la

BUILT_SOURCES = \
	hicum2.c	\
	hicum2.hxx	\
	hicum2acld.c	\
	hicum2ask.c	\
	hicum2defs.h	\
	hicum2del.c	\
	hicum2dest.c	\
	hicum2ext.h	\
	hicum2guesstopology.c	\
	hicum2init.c	\
	hicum2init.h	\
	hicum2itf.h	\
	hicum2load.c	\
	hicum2mask.c	\
	hicum2mdel.c	\
	hicum2mpar.c	\
	hicum2par.c	\
	hicum2pzld.c	\
	hicum2setup.c	\
	hicum2temp.c	\
	hicum2trunc.c

libhicum2_la_SOURCES = 	\
	$(BUILT_SOURCES)

CLEANFILES = \
	$(BUILT_SOURCES)	\
	hicum2noise.c	\
	.hicum2.va.adms \
	.adms.implicit.xml \
	.interface.xml \
	*.h \
	*.xml

#TODO (not implemented) \
	hicum2conv.c	\
	hicum2getic.c

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include
AM_CFLAGS = $(STATIC)

MAINTAINERCLEANFILES = Makefile.in \
	Makefile.am

DISTCLEANFILES = Makefile.am \
	Makefile.in

%.c %.hxx \
	%acld.c %ask.c %defs.h %del.c %dest.c %ext.h %guesstopology.c \
	%init.c %init.h %itf.h %load.c %mask.c %mdel.c %mpar.c %par.c \
	%pzld.c %setup.c %temp.c %trunc.c \
 : $(srcdir)/admsva/%.va \
	$(ADMSXMLINTERFACE)/ngspiceVersion.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEitf.h.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEinit.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEinit.h.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEext.h.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEdefs.h.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEask.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEmask.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEpar.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEmpar.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEload.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEacld.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEpzld.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEtemp.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEtrunc.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEsetup.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEdel.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEmdel.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEdest.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEnoise.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULEguesstopology.c.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULE.hxx.xml \
	$(ADMSXMLINTERFACE)/ngspiceMODULE.c.xml
	admsXml -I$(srcdir) -I$(srcdir)/admsva $< \
	-e $(ADMSXMLINTERFACE)/ngspiceVersion.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEitf.h.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEinit.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEinit.h.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEext.h.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEdefs.h.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEask.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEmask.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEpar.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEmpar.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEload.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEacld.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEpzld.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEtemp.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEtrunc.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEsetup.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEdel.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEmdel.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEdest.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEnoise.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULEguesstopology.c.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULE.hxx.xml \
	-e $(ADMSXMLINTERFACE)/ngspiceMODULE.c.xml
