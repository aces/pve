dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([pve.c])
dnl AC_CONFIG_AUX_DIR(ac_config_aux)
AM_INIT_AUTOMAKE(pve, 0.90)
AM_CONFIG_HEADER([config.h])

AC_REVISION($Revision$)

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CC

smr_WITH_BUILD_PATH

# Check for required C libraries
mni_REQUIRE_VOLUMEIO

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

