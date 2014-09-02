dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([pve.c])
dnl AC_CONFIG_AUX_DIR(ac_config_aux)
AM_INIT_AUTOMAKE(pve, 1.4.2)
AM_CONFIG_HEADER([config.h])

AC_REVISION($Revision$)

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CC
AC_PROG_CXX
smr_WITH_BUILD_PATH

# Check for required C libraries
mni_REQUIRE_VOLUMEIO

# check for perl
AC_PATH_PROGS(PERL, [perl] )

AC_CONFIG_FILES([Makefile scripts/Makefile models/Makefile])
AC_OUTPUT(epm-header
scripts/discretize_pve
scripts/pve_4classes
scripts/pve_curvature
scripts/pve_insect
scripts/pve_script
)

