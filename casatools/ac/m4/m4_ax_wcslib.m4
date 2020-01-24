#
#   AX_WCSLIB([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#

AC_DEFUN([AX_WCSLIB],[

ax_wcslib_ok=no

PKG_CHECK_MODULES(WCSLIB, wcslib, [ax_wcslib_ok=yes], [ax_wcslib_ok=no])

if test x"$ax_wcslib_ok" = x"no"; then
   if test "x$WCSLIB_LIBS" = "x" ; then
      WCSLIB_LIBS="-lwcs"
      AC_ARG_WITH([wcslib-libdir],
	 [  --with-wcslib-libdir=DIR directory where the library was installed],
	 [WCSLIB_LIBS="-L$withval $WCSLIB_LIBS"], )
   fi
   if test "x$WCSLIB_CFLAGS" = "x" ; then
      WCSLIB_CFLAGS=""
      AC_ARG_WITH(wcslib-includedir,
	 [  --with-wcslib-includedir=DIR directory where the headers were installed],
	 [WCSLIB_CFLAGS="-I$withval"], )
   fi
   ax_wcslib_ok=no
   LIBS_sav="$LIBS"
   LIBS="$LIBS $WCSLIB_LIBS"
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $WCSLIB_CFLAGS"
   AC_TRY_COMPILE( [#include <wcslib/wcs.h>],
                   [int v[3]; wcslib_version(v);], 
                   [ax_wcslib_ok=yes],
                   [AC_MSG_WARN([  *** could not find wcslib])] )

   AC_SUBST(WCSLIB_CFLAGS)
   AC_SUBST(WCSLIB_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_wcslib_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_WCSLIB, [1], [Define if you have WCSLIB library.]),[$1])
        :
else
        ax_wcslib_ok=no
        $2
fi
])
