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

AC_DEFUN([AX_OPENMP_CLEAR_CACHE],[dnl
    ## clear openmp here here here here here here here here here
    unset ax_cv_cxx_openmp
])

AC_DEFUN([AX_CXX11_RESET_CHECK],[dnl
    dnl unset cache variables...
    AX_CXX11_CLEAR_CACHE
    AX_OPENMP_CLEAR_CACHE

    ax_found="no"
    ax_new_unchecked_path=""
    ax_cxx11_save_ifs=$IFS
    IFS=$PATH_SEPARATOR
    for ax_dir in $ax_unchecked_path; do
        IFS=$ax_cxx11_save_ifs
        if [[[ "$ax_found" = "yes" ]]]; then
            if [[[ "${ax_new_unchecked_path}" = "" ]]]; then
                ax_new_unchecked_path="${ax_dir}"
            else
                ax_new_unchecked_path="${ax_new_unchecked_path}${PATH_SEPARATOR}${ax_dir}"
            fi
        else
            if [[[ -x "${ax_dir}/${CXX}" ]]]; then
                ax_found="yes"
            fi
            if [[[ "${ax_checked_path}" = "" ]]]; then
                ax_checked_path="${ax_dir}"
            else
                ax_checked_path="${ax_checked_path}${PATH_SEPARATOR}${ax_dir}"
            fi
        fi
    done
    IFS=$ax_cxx11_save_ifs
    ax_unchecked_path=${ax_new_unchecked_path}

])

AC_DEFUN([AX_CXX11], [dnl
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([AX_OPENMP])
    m4_if([$1], [], [m4_fatal([first argument to AX_CXX11 missing])])dnl
    ax_path_saved=$PATH
    ax_checked_path=""
    ax_unchecked_path=$1
    ax_done_cxx11=no
    if [[ "$host_osname" = "darwin" ]]; then
        ax_compiler_choices="clang++ c++ g++ g++-mp-5"
    else
        ax_compiler_choices="g++ c++ clang++"
    fi
    ax_cxx_first_cxx_11=""
    while [[ "$ax_unchecked_path" != "" -a "$ax_done_cxx11" = "no" ]]; do
        PATH=${ax_unchecked_path}${PATH_SEPARATOR}${ax_checked_path}
        for compiler in $ax_compiler_choices; do
            AX_CXX11_CLEAR_CACHE
            AX_OPENMP_CLEAR_CACHE
            ac_save_CXX=$CXX
            CXX="$compiler"
            AC_PROG_CXX($compiler)
            AX_CXX_COMPILE_STDCXX_11(noext,optional)
            if [[ "${ac_success}" = "yes" ]]; then
                if [[ "${ax_cxx_first_cxx_11}" = "" ]]; then
                    ax_cxx_first_cxx_11="$CXX"
                fi
                AX_OPENMP([ax_done_cxx11="yes"; break])
            fi
            CXX="$ac_save_CXX"
        done
        if [[ "$ax_done_cxx11" = "no" ]]; then
            AX_CXX11_RESET_CHECK
        fi
    done
    PATH=$ax_path_saved
    if [[ "${ac_success}" = "yes" ]]; then
        AX_PATH_TO_BINARY(${CXX},CXX,$ax_unchecked_path)
    else
        AC_MSG_WARN([could not find a compiler which supports OpenMP])
        AX_CXX11_RESET_CHECK
        AC_PROG_CXX($ax_cxx_first_cxx_11)
        AC_PATH_PROG([CXX],[$ax_cxx_first_cxx_11])
        ac_success="yes"
    fi
])
