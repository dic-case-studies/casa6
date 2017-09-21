#
# SYNOPSIS
#
#   AX_SPLIT(string,var,[char])
#
# DESCRIPTION
#
#   split up a string (or path) based upon a character
#
# LICENSE
#
#   Copyright (C) 2015 Associated Universities, Inc. Washington DC, USA.
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
AC_DEFUN([AX_SPLIT], [dnl
  m4_if([$1], [], [m4_fatal([first argument to AX_SPLIT missing])])dnl
  m4_if([$2], [], [m4_fatal([second argument to AX_SPLIT missing])])dnl
  m4_if([$3], [], [ax_char=':'], [ax_char=$3])
  ax_save_IFS=$IFS
  IFS=$ax_char
  $2=`echo $1`
  IFS=$ax_save_IFS
  ac_success=yes
])

