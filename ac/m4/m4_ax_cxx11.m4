dnl#
dnl# SYNOPSIS
dnl#
dnl#   AX_CXX11(path)
dnl#
dnl# DESCRIPTION
dnl#
dnl#   find c++ compiler that support
dnl#
dnl# LICENSE
dnl#
dnl#   Copyright (C) 2014 Associated Universities, Inc. Washington DC, USA.
dnl#
dnl#   Copying and distribution of this file, with or without modification, are
dnl#   permitted in any medium without royalty provided the copyright notice
dnl#   and this notice are preserved. This file is offered as-is, without any
dnl#   warranty.

AC_DEFUN([AX_CXX11_CLEAR_CACHE],[dnl
    unset ac_cv_prog_CXX
    unset ax_cv_cxx_compile_cxx11
    unset ac_cv_prog_ac_ct_CXX
    unset ac_cv_cxx_compiler_gnu
    unset ac_cv_prog_cxx_g
    for switch in -std=c++11 -std=c++0x +std=c++11 "-h std=c++11"; do
        cachevar=`$as_echo "ax_cv_cxx_compile_cxx11_$switch" | $as_tr_sh`
        unset $cachevar
    done
])

AC_DEFUN([AX_CXX11], [dnl
    AC_REQUIRE([AC_PROG_CXX])
    m4_if([$1], [], [m4_fatal([first argument to AX_CXX11 missing])])dnl
    ax_path_saved=$PATH
    ax_checked_path=""
    ax_unchecked_path=$1
    while [[ "$ax_unchecked_path" != "" ]]; do
        PATH=${ax_unchecked_path}${PATH_SEPARATOR}${ax_checked_path}
        AC_PROG_CXX(m4_if($host_osname,[darwin],[clang++ c++ g++],[g++ c++ clang++]))
        AX_CXX_COMPILE_STDCXX_11(noext,optional)
        if [[ "${ac_success}" = "yes" ]]; then
            break
        else
            dnl unset cache variables...
            AX_CXX11_CLEAR_CACHE

            ax_found="no"
            ax_new_unchecked_path=""
            ax_cxx11_save_ifs=$IFS
            IFS=$PATH_SEPARATOR
            for ax_dir in $ax_unchecked_path; do
                IFS=$ax_cxx11_save_ifs
                if [[ "$ax_found" = "yes" ]]; then
                    if [[ "${ax_new_unchecked_path}" = "" ]]; then
                        ax_new_unchecked_path="${ax_dir}"
                    else
                        ax_new_unchecked_path="${ax_new_unchecked_path}${PATH_SEPARATOR}${ax_dir}"
                    fi
                else
                    if [[ -x "${ax_dir}/${CXX}" ]]; then
                        ax_found="yes"
                    fi
                    if [[ "${ax_checked_path}" = "" ]]; then
                        ax_checked_path="${ax_dir}"
                    else
                        ax_checked_path="${ax_checked_path}${PATH_SEPARATOR}${ax_dir}"
                    fi
                fi
            done
            IFS=$ax_cxx11_save_ifs
            ax_unchecked_path=${ax_new_unchecked_path}
        fi
    done

    if [[ "${ac_success}" = "yes" ]]; then
        AX_PATH_TO_BINARY(${CXX},CXX,$ax_unchecked_path)
    fi
    PATH=$ax_path_saved
])
