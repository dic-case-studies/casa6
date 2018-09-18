#
# SYNOPSIS
#
#   AX_GSL([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   sets GSL_CFLAGS and GSL_LIBS

AC_DEFUN([AX_GSL],[

  ax_gsl_ok=no
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)

  AC_PATH_PROG([ax_gsl_config], [gsl-config])
  AS_IF([test "x$ax_gsl_config" = "x"],[
      ## unset cache flag
      unset ac_cv_path_ax_gsl_config
      PATH_save="$PATH"
      PATH="$PATH:/opt/casa/03/bin:/opt/casa/02/bin:/opt/local/bin:/usr/local/bin"
      AC_PATH_PROG([ax_gsl_config], [gsl-config])
      PATH="$PATH_save"
  ])

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
