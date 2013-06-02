# Generated automatically from Makefile.in by configure.
# Makefile.in generated automatically by automake 1.4 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

####### kdevelop will overwrite this part!!! (begin)##########


SHELL = /bin/sh

srcdir = .
top_srcdir = ..
prefix = /usr/local
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${prefix}/share
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${prefix}/info
mandir = ${prefix}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/kepler
pkglibdir = $(libdir)/kepler
pkgincludedir = $(includedir)/kepler

top_builddir = ..

ACLOCAL = aclocal
AUTOCONF = autoconf
AUTOMAKE = automake
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL_PROGRAM}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
CC = gcc
CXX = c++
MAKEINFO = makeinfo
PACKAGE = kepler
RANLIB = ranlib
VERSION = 0.1

####### kdevelop will overwrite this part!!! (end)############
bin_PROGRAMS = kepler
kepler_SOURCES = texture.cpp physdat.cpp og_comm.cpp vsop87.cpp tscheby.cpp tp_pos.cpp ppos_lo4.cpp ppos_lo3.cpp ppos_lo2.cpp ppos_lo1.cpp ppos_hi.cpp position.cpp pluto_hi.cpp orbit.cpp monddat.cpp mond_hi.cpp einaus.cpp cartes.cpp astromat.cpp kepler.cpp main.cpp 
kepler_LDADD = -L/usr/X11R6/lib  -lglut -lGLU -lGL -lXmu -lXt -lSM -lICE -lXi  -lm -lX11 -lm -lXext -lX11

SUBDIRS = docs texte texttools 

EXTRA_DIST = main.cpp kepler.h kepler.cpp astromat.cpp cartes.cpp einaus.cpp mond_hi.cpp monddat.cpp orbit.cpp pluto_hi.cpp position.cpp ppos_hi.cpp ppos_lo1.cpp ppos_lo2.cpp ppos_lo3.cpp ppos_lo4.cpp tp_pos.cpp tscheby.cpp vsop87.cpp astromat.h cartes.h definiti.h einaus.h enums.h monddat.h orbit.h position.h ppos_lo.h tp_pos.h tscheby.h vsop87.h vsop87_size.h og_comm.cpp physdat.cpp texture.cpp gv_comm.h physdat.h texture.h erde.ele halley.ele mars.ele merkur.ele planet1.ele planet2.ele venus.ele CropTitel.rgb CropWeiter_DE.rgb CropWeiter_EN.rgb Text0_DE.rgb Text0_EN.rgb Text1_DE.rgb Text1_EN.rgb Text2_DE.rgb Text2_EN.rgb Text3_DE.rgb Text3_EN.rgb Text4_DE.rgb Text4_EN.rgb Text5_DE.rgb Text5_EN.rgb Untertitel1_0_DE.rgb Untertitel1_0_EN.rgb Untertitel1_1_DE.rgb Untertitel1_1_EN.rgb Untertitel1_2_DE.rgb Untertitel1_2_EN.rgb Untertitel2_0_DE.rgb Untertitel2_0_EN.rgb Untertitel3_0_DE.rgb Untertitel3_0_EN.rgb Untertitel3_1_DE.rgb Untertitel3_1_EN.rgb Untertitel3_2_DE.rgb Untertitel3_2_EN.rgb Untertitel4_0_DE.rgb Untertitel4_0_EN.rgb Untertitel5_0_DE.rgb Untertitel5_0_EN.rgb Untertitel5_1_DE.rgb Untertitel5_1_EN.rgb f2 
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_HEADER = ../config.h
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(bin_PROGRAMS)


DEFS = -DHAVE_CONFIG_H -I. -I$(srcdir) -I..
CPPFLAGS = 
LDFLAGS = 
LIBS = 
kepler_OBJECTS =  texture.o physdat.o og_comm.o vsop87.o tscheby.o \
tp_pos.o ppos_lo4.o ppos_lo3.o ppos_lo2.o ppos_lo1.o ppos_hi.o \
position.o pluto_hi.o orbit.o monddat.o mond_hi.o einaus.o cartes.o \
astromat.o kepler.o main.o
kepler_DEPENDENCIES = 
kepler_LDFLAGS = 
CXXFLAGS = -g -O2
CXXCOMPILE = $(CXX) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)
CXXLD = $(CXX)
CXXLINK = $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@
DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(kepler_SOURCES)
OBJECTS = $(kepler_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .cpp .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps kepler/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-binPROGRAMS:

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

distclean-binPROGRAMS:

maintainer-clean-binPROGRAMS:

install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo "  $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`"; \
	     $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	  else :; fi; \
	done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	list='$(bin_PROGRAMS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	done

.c.o:
	$(COMPILE) -c $<

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	-rm -f *.o core *.core

clean-compile:

distclean-compile:
	-rm -f *.tab.c

maintainer-clean-compile:

kepler: $(kepler_OBJECTS) $(kepler_DEPENDENCIES)
	@rm -f kepler
	$(CXXLINK) $(kepler_LDFLAGS) $(kepler_OBJECTS) $(kepler_LDADD) $(LIBS)
.cpp.o:
	$(CXXCOMPILE) -c $<

# This directory's subdirectories are mostly independent; you can cd
# into them and run `make' without going through this Makefile.
# To change the values of `make' variables: instead of editing Makefiles,
# (1) if the variable is set in `config.status', edit `config.status'
#     (which will cause the Makefiles to be regenerated when you run `make');
# (2) otherwise, pass the desired values on the `make' command line.



all-recursive install-data-recursive install-exec-recursive \
installdirs-recursive install-recursive uninstall-recursive  \
check-recursive installcheck-recursive info-recursive dvi-recursive:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	dot_seen=no; \
	target=`echo $@ | sed s/-recursive//`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    dot_seen=yes; \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	   || case "$$amf" in *=*) exit 1;; *k*) fail=yes;; *) exit 1;; esac; \
	done; \
	if test "$$dot_seen" = "no"; then \
	  $(MAKE) $(AM_MAKEFLAGS) "$$target-am" || exit 1; \
	fi; test -z "$$fail"

mostlyclean-recursive clean-recursive distclean-recursive \
maintainer-clean-recursive:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	dot_seen=no; \
	rev=''; list='$(SUBDIRS)'; for subdir in $$list; do \
	  rev="$$subdir $$rev"; \
	  test "$$subdir" = "." && dot_seen=yes; \
	done; \
	test "$$dot_seen" = "no" && rev=". $$rev"; \
	target=`echo $@ | sed s/-recursive//`; \
	for subdir in $$rev; do \
	  echo "Making $$target in $$subdir"; \
	  if test "$$subdir" = "."; then \
	    local_target="$$target-am"; \
	  else \
	    local_target="$$target"; \
	  fi; \
	  (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) $$local_target) \
	   || case "$$amf" in *=*) exit 1;; *k*) fail=yes;; *) exit 1;; esac; \
	done && test -z "$$fail"
tags-recursive:
	list='$(SUBDIRS)'; for subdir in $$list; do \
	  test "$$subdir" = . || (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) tags); \
	done

tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	here=`pwd` && cd $(srcdir) \
	  && mkid -f$$here/ID $$unique $(LISP)

TAGS: tags-recursive $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
   if test "$$subdir" = .; then :; else \
	    test -f $$subdir/TAGS && tags="$$tags -i $$here/$$subdir/TAGS"; \
   fi; \
	done; \
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)$$unique$(LISP)$$tags" \
	  || (cd $(srcdir) && etags $(ETAGS_ARGS) $$tags  $$unique $(LISP) -o $$here/TAGS)

mostlyclean-tags:

clean-tags:

distclean-tags:
	-rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = kepler

distdir: $(DISTFILES)
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$d/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done
	for subdir in $(SUBDIRS); do \
	  if test "$$subdir" = .; then :; else \
	    test -d $(distdir)/$$subdir \
	    || mkdir $(distdir)/$$subdir \
	    || exit 1; \
	    chmod 777 $(distdir)/$$subdir; \
	    (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) top_distdir=../$(top_distdir) distdir=../$(distdir)/$$subdir distdir) \
	      || exit 1; \
	  fi; \
	done
astromat.o: astromat.cpp astromat.h cartes.h
cartes.o: cartes.cpp cartes.h
einaus.o: einaus.cpp astromat.h cartes.h einaus.h enums.h
kepler.o: kepler.cpp astromat.h cartes.h definiti.h enums.h physdat.h \
	orbit.h position.h ppos_lo.h tscheby.h gv_comm.h kepler.h
main.o: main.cpp ../config.h astromat.h cartes.h einaus.h enums.h \
	definiti.h physdat.h orbit.h position.h ppos_lo.h tscheby.h \
	monddat.h gv_comm.h texture.h kepler.h
mond_hi.o: mond_hi.cpp astromat.h cartes.h
monddat.o: monddat.cpp astromat.h cartes.h monddat.h
og_comm.o: og_comm.cpp astromat.h cartes.h definiti.h enums.h physdat.h \
	orbit.h position.h ppos_lo.h tscheby.h gv_comm.h
orbit.o: orbit.cpp astromat.h cartes.h orbit.h
physdat.o: physdat.cpp physdat.h astromat.h cartes.h enums.h
pluto_hi.o: pluto_hi.cpp astromat.h cartes.h einaus.h enums.h
position.o: position.cpp astromat.h cartes.h orbit.h definiti.h enums.h \
	position.h ppos_lo.h
ppos_hi.o: ppos_hi.cpp astromat.h cartes.h definiti.h enums.h vsop87.h
ppos_lo1.o: ppos_lo1.cpp ppos_lo.h astromat.h cartes.h
ppos_lo2.o: ppos_lo2.cpp ppos_lo.h astromat.h cartes.h
ppos_lo3.o: ppos_lo3.cpp ppos_lo.h astromat.h cartes.h
ppos_lo4.o: ppos_lo4.cpp ppos_lo.h astromat.h cartes.h
texture.o: texture.cpp
tp_pos.o: tp_pos.cpp astromat.h cartes.h definiti.h enums.h position.h \
	ppos_lo.h tscheby.h tp_pos.h
tscheby.o: tscheby.cpp astromat.h cartes.h definiti.h enums.h tscheby.h
vsop87.o: vsop87.cpp astromat.h cartes.h vsop87.h vsop87_size.h

info-am:
info: info-recursive
dvi-am:
dvi: dvi-recursive
check-am: all-am
check: check-recursive
installcheck-am:
installcheck: installcheck-recursive
install-exec-am: install-binPROGRAMS
install-exec: install-exec-recursive

install-data-am: install-data-local
install-data: install-data-recursive

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-recursive
uninstall-am: uninstall-binPROGRAMS uninstall-local
uninstall: uninstall-recursive
all-am: Makefile $(PROGRAMS)
all-redirect: all-recursive
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs: installdirs-recursive
installdirs-am:
	$(mkinstalldirs)  $(DESTDIR)$(bindir)


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
mostlyclean-am:  mostlyclean-binPROGRAMS mostlyclean-compile \
		mostlyclean-tags mostlyclean-generic

mostlyclean: mostlyclean-recursive

clean-am:  clean-binPROGRAMS clean-compile clean-tags clean-generic \
		mostlyclean-am

clean: clean-recursive

distclean-am:  distclean-binPROGRAMS distclean-compile distclean-tags \
		distclean-generic clean-am

distclean: distclean-recursive

maintainer-clean-am:  maintainer-clean-binPROGRAMS \
		maintainer-clean-compile maintainer-clean-tags \
		maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-recursive

.PHONY: mostlyclean-binPROGRAMS distclean-binPROGRAMS clean-binPROGRAMS \
maintainer-clean-binPROGRAMS uninstall-binPROGRAMS install-binPROGRAMS \
mostlyclean-compile distclean-compile clean-compile \
maintainer-clean-compile install-data-recursive \
uninstall-data-recursive install-exec-recursive \
uninstall-exec-recursive installdirs-recursive uninstalldirs-recursive \
all-recursive check-recursive installcheck-recursive info-recursive \
dvi-recursive mostlyclean-recursive distclean-recursive clean-recursive \
maintainer-clean-recursive tags tags-recursive mostlyclean-tags \
distclean-tags clean-tags maintainer-clean-tags distdir info-am info \
dvi-am dvi check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-local install-data-am install-data install-am \
install uninstall-local uninstall-am uninstall all-redirect all-am all \
installdirs-am installdirs mostlyclean-generic distclean-generic \
clean-generic maintainer-clean-generic clean mostlyclean distclean \
maintainer-clean


install-data-local:
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) erde.ele ${bindir}/erde.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) halley.ele ${bindir}/halley.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) mars.ele ${bindir}/mars.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) merkur.ele ${bindir}/merkur.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) planet1.ele ${bindir}/planet1.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) planet2.ele ${bindir}/planet2.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) venus.ele ${bindir}/venus.ele
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) CropTitel.rgb ${bindir}/CropTitel.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) CropWeiter_DE.rgb ${bindir}/CropWeiter_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) CropWeiter_EN.rgb ${bindir}/CropWeiter_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text0_DE.rgb ${bindir}/Text0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text0_EN.rgb ${bindir}/Text0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text1_DE.rgb ${bindir}/Text1_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text1_EN.rgb ${bindir}/Text1_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text2_DE.rgb ${bindir}/Text2_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text2_EN.rgb ${bindir}/Text2_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text3_DE.rgb ${bindir}/Text3_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text3_EN.rgb ${bindir}/Text3_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text4_DE.rgb ${bindir}/Text4_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text4_EN.rgb ${bindir}/Text4_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text5_DE.rgb ${bindir}/Text5_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Text5_EN.rgb ${bindir}/Text5_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_0_DE.rgb ${bindir}/Untertitel1_0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_0_EN.rgb ${bindir}/Untertitel1_0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_1_DE.rgb ${bindir}/Untertitel1_1_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_1_EN.rgb ${bindir}/Untertitel1_1_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_2_DE.rgb ${bindir}/Untertitel1_2_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel1_2_EN.rgb ${bindir}/Untertitel1_2_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel2_0_DE.rgb ${bindir}/Untertitel2_0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel2_0_EN.rgb ${bindir}/Untertitel2_0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_0_DE.rgb ${bindir}/Untertitel3_0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_0_EN.rgb ${bindir}/Untertitel3_0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_1_DE.rgb ${bindir}/Untertitel3_1_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_1_EN.rgb ${bindir}/Untertitel3_1_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_2_DE.rgb ${bindir}/Untertitel3_2_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel3_2_EN.rgb ${bindir}/Untertitel3_2_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel4_0_DE.rgb ${bindir}/Untertitel4_0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel4_0_EN.rgb ${bindir}/Untertitel4_0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel5_0_DE.rgb ${bindir}/Untertitel5_0_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel5_0_EN.rgb ${bindir}/Untertitel5_0_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel5_1_DE.rgb ${bindir}/Untertitel5_1_DE.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) Untertitel5_1_EN.rgb ${bindir}/Untertitel5_1_EN.rgb
	$(mkinstalldirs) ${bindir}/
	$(INSTALL_DATA) f2 ${bindir}/f2

uninstall-local:
	-rm -f ${bindir}/erde.ele
	-rm -f ${bindir}/halley.ele
	-rm -f ${bindir}/mars.ele
	-rm -f ${bindir}/merkur.ele
	-rm -f ${bindir}/planet1.ele
	-rm -f ${bindir}/planet2.ele
	-rm -f ${bindir}/venus.ele
	-rm -f ${bindir}/CropTitel.rgb
	-rm -f ${bindir}/CropWeiter_DE.rgb
	-rm -f ${bindir}/CropWeiter_EN.rgb
	-rm -f ${bindir}/Text0_DE.rgb
	-rm -f ${bindir}/Text0_EN.rgb
	-rm -f ${bindir}/Text1_DE.rgb
	-rm -f ${bindir}/Text1_EN.rgb
	-rm -f ${bindir}/Text2_DE.rgb
	-rm -f ${bindir}/Text2_EN.rgb
	-rm -f ${bindir}/Text3_DE.rgb
	-rm -f ${bindir}/Text3_EN.rgb
	-rm -f ${bindir}/Text4_DE.rgb
	-rm -f ${bindir}/Text4_EN.rgb
	-rm -f ${bindir}/Text5_DE.rgb
	-rm -f ${bindir}/Text5_EN.rgb
	-rm -f ${bindir}/Untertitel1_0_DE.rgb
	-rm -f ${bindir}/Untertitel1_0_EN.rgb
	-rm -f ${bindir}/Untertitel1_1_DE.rgb
	-rm -f ${bindir}/Untertitel1_1_EN.rgb
	-rm -f ${bindir}/Untertitel1_2_DE.rgb
	-rm -f ${bindir}/Untertitel1_2_EN.rgb
	-rm -f ${bindir}/Untertitel2_0_DE.rgb
	-rm -f ${bindir}/Untertitel2_0_EN.rgb
	-rm -f ${bindir}/Untertitel3_0_DE.rgb
	-rm -f ${bindir}/Untertitel3_0_EN.rgb
	-rm -f ${bindir}/Untertitel3_1_DE.rgb
	-rm -f ${bindir}/Untertitel3_1_EN.rgb
	-rm -f ${bindir}/Untertitel3_2_DE.rgb
	-rm -f ${bindir}/Untertitel3_2_EN.rgb
	-rm -f ${bindir}/Untertitel4_0_DE.rgb
	-rm -f ${bindir}/Untertitel4_0_EN.rgb
	-rm -f ${bindir}/Untertitel5_0_DE.rgb
	-rm -f ${bindir}/Untertitel5_0_EN.rgb
	-rm -f ${bindir}/Untertitel5_1_DE.rgb
	-rm -f ${bindir}/Untertitel5_1_EN.rgb
	-rm -f ${bindir}/f2

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
