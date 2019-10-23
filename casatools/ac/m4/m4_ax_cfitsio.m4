#
# SYNOPSIS
#
#   AX_CFITSIO([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for the cfitsio library (see
#   http://heasarc.gsfc.nasa.gov/fitsio/). On success, it sets the
#   CFITSIO_LIBS and CFITSIO_CFLAGS output variables to hold the
#   requisite library linkages. It first look for a pkg-config file
#   before trying to compile and link with cfitsio.
#
#   To link with CFITSIO, you should link with:
#
#     $CFITSIO_LIBS $LIBS
#
#   ACTION-IF-FOUND is a list of shell commands to run if the CFITSIO
#   library is found, and ACTION-IF-NOT-FOUND is a list of commands to
#   run it if it is not found. If ACTION-IF-FOUND is not specified, the
#   default action will define the processor macro HAVE_CFITSIO. 
#
# LICENSE
#
#   Copyright (c) 2010 Sebastien Fabbro <sfabbro@uvic.ca>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.


AC_DEFUN([AX_CFITSIO],[

ax_cfitsio_ok=no

PKG_CHECK_MODULES(CFITSIO, cfitsio, [ax_cfitsio_ok=yes], [ax_cfitsio_ok=no])

###
### if PKG_CHECK_MODULES fails, this code is of little help...
### it seems like it needs to:
###        (1) find the default libcfitsio.so
###        (2) ascertain it's fully qualified path
###        (3) use the path to the library as the starting
###            point for finding fitsio.h
###
if test x"$ax_cfitsio_ok" = x"no"; then
   # if neither pkg-config file found, nor CFITSIO_LIBS defined
   if test x$CFITSIO_LIBS = x ; then
      CFITSIO_LIBS="-lcfitsio -lm"
      AC_ARG_WITH([cfitsio-libdir],
	 [  --with-cfitsio-libdir=DIR directory where the library was installed],
	 [CFITSIO_LIBS="-L$withval $CFITSIO_LIBS"], )
   fi
   ax_cfitsio_lib_ok=no
   LIBS_sav="$LIBS"
   LIBS="$LIBS $CFITSIO_LIBS"
   AC_CHECK_LIB([cfitsio], [ffopen], [ax_cfitsio_lib_ok=yes], [AC_MSG_WARN([  *** cfitsio library not found])])

   # if neither pkg-config file found, nor CFITSIO_CFLAGS defined
   if test x$CFITSIO_CFLAGS = x ; then
      CFITSIO_CFLAGS=""
      AC_ARG_WITH(cfitsio-includedir,
	 [  --with-cfitsio-includedir=DIR directory where the headers were installed],
	 [CFITSIO_CFLAGS="-I$withval"], )
   fi
   ax_cfitsio_hdr_ok=no
   CPPFLAGS_sav="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $CFITSIO_CFLAGS"
   AC_CHECK_HEADER([fitsio.h],
      [ax_cfitsio_hdr_ok=yes],
      [AC_MSG_WARN([  *** cfitsio headers not found.])])

   if test x$ax_cfitsio_lib_ok = xyes -a x$ax_cfitsio_hdr_ok = xyes; then
      ax_cfitsio_ok=yes
   fi
   AC_SUBST(CFITSIO_CFLAGS)
   AC_SUBST(CFITSIO_LIBS)
   LIBS="$LIBS_sav"
   CPPFLAGS="$CPPFLAGS_sav"
fi


# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_cfitsio_ok" = x"yes"; then
        ifelse([$1],,AC_DEFINE(HAVE_CFITSIO, [1], [Define if you have CFITSIO library.]),[$1])
        :
else
        ax_cfitsio_ok=no
        $2
fi
])
