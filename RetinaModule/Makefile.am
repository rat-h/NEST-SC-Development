# Automake file for external dynamic modules for NEST
#
# Hans Ekkehard Plesser, April 2008
# Automake file for the Developer Module
#
# libretinamodule is built as a normal, installable library.
# It will be installed to $prefix/lib by make install.
#
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.


# 1. Exchange "my" in "retinamodule" with the name of your model below
#    (ten locations).
# 2. Add all .cpp and .h files from your code as *_SOURCES. Header files
#    are given only so that they will be included in the tarball if you
#    run "make dist" on your module.
# 3. The retinamodule* stuff creates a module that can be loaded at runtime.
#    It is called retinamodule.so.
# 4. The libretinamodule* stuff creates a library against which NEST can be
#    linked.

libdir= @libdir@/nest

# We need to install the module header for static linking on BlueGene
include_HEADERS = retinamodule.h

# All other source files
source_files=  retinamodule.cpp \
               stbrst_gc_conv.cpp stbrst_gc_conv.h \
               convolution.h

if BUILD_DYNAMIC_USER_MODULES
  lib_LTLIBRARIES= libretinamodule.la retinamodule.la
  retinamodule_la_CXXFLAGS= @AM_CXXFLAGS@
  retinamodule_la_SOURCES=  $(source_files)
  retinamodule_la_LDFLAGS=  -module
else
  lib_LTLIBRARIES= libretinamodule.la
endif

libretinamodule_la_CXXFLAGS= @AM_CXXFLAGS@ -DLINKED_MODULE
libretinamodule_la_SOURCES=  $(source_files)

MAKEFLAGS= @MAKE_FLAGS@

AM_CPPFLAGS= @NEST_CPPFLAGS@ \
             @INCLTDL@

.PHONY: install-slidoc

pkgdatadir=@datadir@/nest

nobase_pkgdata_DATA=\
	sli/retinamodule-init.sli

install-slidoc:
	NESTRCFILENAME=/dev/null $(DESTDIR)$(NEST_PREFIX)/bin/sli --userargs="@HELPDIRS@" $(NEST_PREFIX)/share/nest/sli/install-help.sli

install-data-hook: install-exec install-slidoc

EXTRA_DIST= sli
