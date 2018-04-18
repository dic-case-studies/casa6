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
    AC_MSG_ERROR([no gsl-config executable found]) 
  ],[
    ax_gsl_uname=`uname`
    AS_IF([test "x$ax_gsl_uname" = "xLinux"],[
      ax_gsl_path=`${ax_gsl_config} --prefix`
      AS_IF([test -e "$ax_gsl_path/lib64/libgsl.so"],[
        ax_gsl_lflag="-L$ax_gsl_path/lib64"
      ],[
        AS_IF([test -e "$ax_gsl_path/lib/libgsl.so"],[
          ax_gsl_lflag="-L$ax_gsl_path/lib"
        ],[ax_gsl_lflag=""])
      ])
    ])
    GSL_LIBS=`${ax_gsl_config} --libs`
    AS_IF([test "x$ax_gsl_lflag" != "x"],[
      GSL_LIBS="${ax_gsl_lflag} ${GSL_LIBS}"
    ])
    GSL_CFLAGS=`${ax_gsl_config} --cflags`
  ])
])
