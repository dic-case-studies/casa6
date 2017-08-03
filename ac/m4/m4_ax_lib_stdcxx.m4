#
# SYNOPSIS
#
#   AX_LIB_STDCXX
#
# DESCRIPTION
#
#   Check for flags needed for new libc++ support.
#
#   This is somewhat dependent on the internals of the generated configure file.
#   $ac_link is manipulated to decouple the compiling from the linking. This
#   allows AC_TRY_COMPILE to be used to generate the object file thereby checking
#   for the requirement of particular flags when compiling, i.e. CXXFLAGS, while
#   AC_LINK_IFELSE is used (via a munged $ac_link) to test the requirement of
#   a flag at link time, i.e. LDFLAGS (e.g. when linking a number of object files
#   without C++ source code).
#
#   Aside: probably CPPFLAGS is the right place to put compile time flags because
#          probably selection of the C++ library entails different sets of header
#          files, but for some compilers, CXXFLAGS may be the right place
#
# LICENSE
#
#   Copyright (C) 2014 Associated Universities, Inc. Washington DC, USA.
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
#
AC_DEFUN([AX_LIB_STDCXX], [
    AC_LANG_PUSH([C++])
    AC_CACHE_CHECK( [for flags to support C++11 libc++], ax_cv_lib_stdcxx_cppflag, [
        AC_CACHE_VAL( ax_cv_lib_stdcxx_ldflag, [
            orig_CPPFLAGS="${CPPFLAGS}"
            orig_LDFLAGS="${LD_FLAGS}"
            ax_cv_lib_stdcxx_cppflag="no"
            ax_cv_lib_stdcxx_ldflag="no"
            for cppflag in "" "-stdlib=libc++"; do
                CPPFLAGS="${orig_CPPFLAGS} ${cppflag}"
                AC_TRY_COMPILE( [#include <memory>],
                                [auto val = std::make_shared<int>(90);],
                                [
                                  ax_cv_lib_stdcxx_cppflag="${cppflag}"
                                  CPPFLAGS="${orig_CPPFLAGS}"
                                  mv conftest.$ac_objext xconftest.$ac_objext
                                  orig_ac_link="$ac_link"
                                  ac_link=`echo $ac_link | sed 's|conftest.$ac_ext|xconftest.$ac_objext|'`
                                  for ldflag in "" "-stdlib=libc++"; do
                                      LDFLAGS="${orig_LDFLAGS} ${ldflag}"
                                      AC_LINK_IFELSE(, [ax_cv_lib_stdcxx_ldflag="${ldflag}"],
                                                       [ax_cv_lib_stdcxx_ldflag=no])
                                  done
                                  ac_link="${orig_ac_link}"
                                  rm -f  xconftest.$ac_objext
                                ],
                                [ax_cv_lib_stdcxx_cppflag=no] )
                if test "${ax_cv_lib_stdcxx_cppflag}" != "no" -a \
                        "${ax_cv_lib_stdcxx_ldflag}" != "no" ; then
                    break
                else
                    ax_cv_lib_stdcxx_cppflag=no
                    ax_cv_lib_stdcxx_ldflag=no
                fi
            done
            if test "${ax_cv_lib_stdcxx_cppflag}" != "no"; then
                CPPFLAGS="${orig_CPPFLAGS} ${ax_cv_lib_stdcxx_cppflag}"
                LDFLAGS="${orig_LDFLAGS} ${ax_cv_lib_stdcxx_ldflag}"
            else
                CPPFLAGS="${orig_CPPFLAGS}"
                LDFLAGS="${orig_LDFLAGS}"
            fi
        ])
    ])
    AC_LANG_POP([C++])
])
