# Makefile.am for MINC
#
AUTOMAKE_OPTIONS = 
ACLOCAL_AMFLAGS = -I m4

# Header files which must be included in the distribution but which are 
# not part of the installation.
#
noinst_HEADERS = \
	least_trimmed_squares.h \
	matrix3.h \
	minvarellipsoid.h \
	minvarellipsoid3.h \
	pve3_aux.h \
	pve_aux.h

# Programs we build and install.
#
bin_PROGRAMS = \
	pve \
	pve3

EXTRA_DIST = $(m4_files)

m4_files = \
	m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4

pve_SOURCES = \
	pve.c \
	pve_aux.c \
	least_trimmed_squares.c \
	minvarellipsoid.c \
	kth_smallest.c

pve3_SOURCES = \
	pve3.c \
	pve3_aux.c \
	matrix3.c \
	minvarellipsoid3.c \
	kth_smallest.c

