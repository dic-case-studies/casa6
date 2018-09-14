#
#   AX_RPFITS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#

AC_DEFUN([AX_RPFITS],[

ax_rpfits_ok=no

if test x"$ax_rpfits_ok" = x"no"; then
   if test "x$RPFITS_LIBS" = "x" ; then
      RPFITS_LIBS="-lrpfits"
      AC_ARG_WITH([rpfits-libdir],
	 [  --with-rpfits-libdir=DIR directory where the library was installed],
	 [RPFITS_LIBS="-L$withval $RPFITS_LIBS"], )
   fi
   if test "x$RPFITS_CFLAGS" = "x" ; then
      RPFITS_CFLAGS=""
      AC_ARG_WITH(rpfits-includedir,
	 [  --with-rpfits-includedir=DIR directory where the headers were installed],
	 [RPFITS_CFLAGS="-I$withval"], )
   fi
   ax_rpfits_ok=no
   LIBS_sav="$LIBS"
   LIBS_sav="$LIBS"
   LIBS="$LIBS $RPFITS_LIBS"
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $RPFITS_CFLAGS"

   AC_CHECK_LIB([rpfits], [rpfitsin_], [ax_rpfits_ok=yes])

   if test "x$ax_rpfits_ok" = "xno"; then
      # fish around...
      for path in [/opt/casa/03/lib /opt/casa/02/lib /opt/local/lib /usr/local/lib]; do
         RPFITS_LIBS="-L$path -lrpfits"
         LIBS="$LIBS_sav $RPFITS_LIBS"
         # prevent the cache from thwarting our efforts...
         unset ac_cv_lib_rpfits_rpfitsin_
         AC_CHECK_LIB([rpfits], [rpfitsin_], [ax_rpfits_ok=yes])
         if test "x$ax_rpfits_ok" = "xyes"; then
            break
         fi
      done
   fi

   if test "x$ax_rpfits_ok" = "xno"; then
      RPFITS_LIBS=""
      RPFITS_CFLAGS=""
   fi

   AC_SUBST(RPFITS_CFLAGS)
   AC_SUBST(RPFITS_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_rpfits_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_RPFITS, [1], [Define if you have RPFITS library.]),[$1])
        :
else
        ax_rpfits_ok=no
        $2
fi
])
