dnl#
dnl# SYNOPSIS
dnl#
dnl#   AX_EXPAND_PATH(path,var)
dnl#
dnl# DESCRIPTION
dnl#
dnl#   expand the path to a fully qualified path and assign it to var
dnl#
dnl# LICENSE
dnl#
dnl#   Copyright (C) 2014,2015 Associated Universities, Inc. Washington DC, USA.
dnl#
dnl#   Copying and distribution of this file, with or without modification, are
dnl#   permitted in any medium without royalty provided the copyright notice
dnl#   and this notice are preserved. This file is offered as-is, without any
dnl#   warranty.

AC_DEFUN([AX_EXPAND_PATH], [dnl
  AC_PREREQ([2.50])dnl
  m4_if([$1], [], [m4_fatal([first argument of AX_EXPAND_PATH missing])])dnl
  m4_if([$2], [],[m4_fatal([second argument of AX_EXPAND_PATH missing])])dnl
  ac_success=no
  if test -d $1; then
    $2=`cd $1 && pwd`
    ac_success=yes
  else
    $2=`cd $(dirname $1) && pwd`/$(basename $1)
    ac_success=yes
  fi
])
