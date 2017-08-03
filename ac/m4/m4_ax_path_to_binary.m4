dnl#
dnl# SYNOPSIS
dnl#
dnl#   AX_PATH_TO_BINARY(executable,var,[PATH])
dnl#
dnl# DESCRIPTION
dnl#
dnl#   find the qualified path to a binary based upon PATH...
dnl#   e.g. expand "gcc" to "/usr/bin/gcc"
dnl#
dnl# LICENSE
dnl#
dnl#   Copyright (C) 2015 Associated Universities, Inc. Washington DC, USA.
dnl#
dnl#   Copying and distribution of this file, with or without modification, are
dnl#   permitted in any medium without royalty provided the copyright notice
dnl#   and this notice are preserved. This file is offered as-is, without any
dnl#   warranty.

AC_DEFUN([AX_PATH_TO_BINARY], [dnl
    m4_if([$1], [], [m4_fatal([first argument of AX_PATH_TO_BINARY missing])])dnl
    m4_if([$2], [],[m4_fatal([second argument of AX_PATH_TO_BINARY missing])])dnl
    ax_path_to_binary_save_ifs=$IFS
    IFS=$PATH_SEPARATOR
    for ax_dir in m4_if([$3],[],[$PATH],[$3]); do
        IFS=$ax_path_to_binary_save_ifs
        if [[ -x "${ax_dir}/$1" ]]; then
            AX_EXPAND_PATH([${ax_dir}/$1],$2)
            break
        fi
    done
    IFS=$ax_path_to_binary_save_ifs
    AC_SUBST($2)
])
