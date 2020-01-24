#
# SYNOPSIS
#
#   AX_GSL([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   sets GSL_CFLAGS and GSL_LIBS

AC_DEFUN([AX_GSL],[

  GSL_REQUEST_VERSION=

  changequote(<<, >>)
  for a in $1 $2 $3 $4 $5 $6 $7 $8 $9 x; do
    case "$a" in
        x) break;;
        [0-9]*.[0-9]*) GSL_REQUEST_VERSION="$a";;
    esac
  done
  changequote([, ])

  changequote(<<, >>)
  gsl_major_req=`expr $GSL_REQUEST_VERSION : '\([0-9]*\)\.[0-9]*.*'`
  gsl_minor_req=`expr $GSL_REQUEST_VERSION : '[0-9]*\.\([0-9]*\).*'`
  changequote([, ])

  ax_gsl_ok=no
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)

  PATH_save="$PATH"
  PATH="$PATH:/opt/casa/03/bin:/opt/casa/02/bin:/opt/local/bin:/usr/local/bin"

  AC_PATH_PROG([ax_gsl_config], [gsl-config])

  AS_IF([test "x$ax_gsl_config" != "x"],[

      AC_MSG_CHECKING(requested gsl version ($GSL_REQUEST_VERSION))

      changequote(<<, >>)
      gsl_version=`${ax_gsl_config} --version 2>&1 | sed 's/.*\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p; d'`
      gsl_major_ver=`echo $gsl_version | sed 's/^\([0-9][0-9]*\).*/\1/'`
      gsl_minor_ver=`echo $gsl_version | sed 's/^[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'`
      changequote([, ])

      if test $gsl_major_ver -ge $gsl_major_req &&
         test $gsl_minor_ver -ge $gsl_minor_req
      then
          AC_MSG_RESULT(yes)
      else
         AC_MSG_RESULT(no)
         while test "x$ax_gsl_config" != "x"; do
             ## save paths to things we need before we ruin PATH
             gsl_sed=`which sed`
             gsl_expr=`which expr`
             gsl_dirname=`which dirname`
             ## unset cache flag
             unset ac_cv_path_ax_gsl_config
             ## save unsuccessful path
             baddir=`$gsl_dirname $ax_gsl_config`
             ## unset gsl config variable
             unset ax_gsl_config
             ## remove unsuccessful path
             old_path=`echo $PATH | $gsl_sed 's|:| |g'`
             PATH=""
             for p in $old_path; do
                 if test "$baddir" != "$p"; then
                     if test "x$PATH" = "x"; then
                         PATH="$p"
                     else
                         PATH="$PATH:$p"
                     fi
                 fi
             done
             AC_PATH_PROG([ax_gsl_config], [gsl-config])   
             if test -n "$ax_gsl_config"; then
                 AC_MSG_CHECKING(requested gsl version ($GSL_REQUEST_VERSION))
                 changequote(<<, >>)
                 gsl_version=`${ax_gsl_config} --version 2>&1 | $gsl_sed 's/.*\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p; d'`
                 gsl_major_ver=`echo $gsl_version | $gsl_sed 's/^\([0-9][0-9]*\).*/\1/'`
                 gsl_minor_ver=`echo $gsl_version | $gsl_sed 's/^[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'`
                 changequote([, ])

                 if test $gsl_major_ver -ge $gsl_major_req &&
                    test $gsl_minor_ver -ge $gsl_minor_req
                 then
                    AC_MSG_RESULT(yes)
                    break
                 else
                    AC_MSG_RESULT(no)
                 fi
             fi
          done
      fi
  ])

  PATH="$PATH_save"

  AS_IF([test "x$ax_gsl_config" = "x"],[
    AC_MSG_ERROR([no gsl-config executable found]) 
  ],[
    ax_gsl_uname=`uname`
    AS_IF([test "x$ax_gsl_uname" = "xLinux"],[
      ax_gsl_path=`${ax_gsl_config} --prefix`
      AS_IF([test -e "$ax_gsl_path/lib64/libgsl.so"],[
        ax_gsl_lflag="-L$ax_gsl_path/lib64"
          ax_gsl_lpath=`${ax_gsl_config} --libs | sed "s|-lgsl\([[a-zA-Z0-9_-]]*\)|${ax_gsl_path}/lib64/libgsl\1.so|g"`
dnl       without an rpath, the GSL library may not be found because it may be 
dnl       outside of the regular dl search path...
          ax_gsl_lpath="${ax_gsl_lpath} -Wl,-rpath,${ax_gsl_path}/lib64"
      ],[
        AS_IF([test -e "$ax_gsl_path/lib/libgsl.so"],[
          ax_gsl_lflag="-L$ax_gsl_path/lib"
          ax_gsl_lpath=`${ax_gsl_config} --libs | sed "s|-lgsl\([[a-zA-Z0-9_-]]*\)|${ax_gsl_path}/lib/libgsl\1.so|g"`
dnl       without an rpath, the GSL library may not be found because it may be 
dnl       outside of the regular dl search path...
          ax_gsl_lpath="${ax_gsl_lpath} -Wl,-rpath,${ax_gsl_path}/lib"
        ],[ax_gsl_lflag="";ax_gsl_lpath=""])
      ])
    ])
    GSL_LIBS=`${ax_gsl_config} --libs`
    AS_IF([test "x$ax_gsl_lpath" != "x"],[
dnl     RedHat Enterprise Linux 7 still has gsl v1.5 libraries in /usr/lib64
dnl     so we go to the extra hassle here of including fully qualified paths
dnl     to the GSL libraries (which for RHEL7 will be somewhere other than
dnl     /usr/lib64 (unfortunately)...
dnl     ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
dnl     GSL_LIBS="${ax_gsl_lflag} ${GSL_LIBS}"
        GSL_LIBS="${ax_gsl_lpath}"
    ])
    GSL_CFLAGS=`${ax_gsl_config} --cflags`
  ])
])
