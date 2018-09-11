dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_pkg_swig.html
dnl
dnl renamed from ac_pkg_swig to ax_swig
AC_DEFUN([AX_SWIG],
[
SWIG_REQUEST_VERSION=

changequote(<<, >>)

for a in $1 $2 $3 $4 $5 $6 $7 $8 $9 x; do
    case "$a" in
        x) break;;
        [0-9]*.[0-9]*.[0-9]*) SWIG_REQUEST_VERSION="$a";;
		c++) SWIGFLAGS="$SWIGFLAGS -c++";;
		raw) SWIGFLAGS="$SWIGFLAGS -c";;
    esac
done

changequote([, ])

AC_PATH_PROG(SWIG,swig)

if test -n "$SWIG";
then
	SWIGLIB=`$SWIG -swiglib`

	AC_SUBST(SWIG)
	AC_SUBST(SWIGLIB)
	AC_SUBST(SWIGFLAGS)

	AC_MSG_CHECKING(swig version)

	changequote(<<, >>)
	swig_version=`$SWIG -version 2>&1 | sed 's/.* \([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p; d'`
	swig_major_ver=`echo $swig_version | sed 's/^\([0-9][0-9]*\).*/\1/'`
	swig_minor_ver=`echo $swig_version | sed 's/^[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'`
	swig_micro_ver=`echo $swig_version | sed 's/^[0-9][0-9]*\.[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'`
	changequote([, ])

	AC_MSG_RESULT($swig_version)

	SWIGVERNUM=`printf "%02d%02d%02d" $swig_major_ver $swig_minor_ver $swig_micro_ver`
	# SWIGVERNUM=`echo $SWIG_REQUEST_VERSION | awk '{ split($[1],a,"\."); print [a[1]*1000000+a[2]*1000+a[3]] }' 2>/dev/null`

	if test -n "$SWIG_REQUEST_VERSION";
	then
		AC_MSG_CHECKING(requested swig version ($SWIG_REQUEST_VERSION))

		changequote(<<, >>)
		swig_major_req=`expr $SWIG_REQUEST_VERSION : '\([0-9]*\)\.[0-9]*\.[0-9]*'`
		swig_minor_req=`expr $SWIG_REQUEST_VERSION : '[0-9]*\.\([0-9]*\)\.[0-9]*'`
		swig_micro_req=`expr $SWIG_REQUEST_VERSION : '[0-9]*\.[0-9]*\.\([0-9]*\)'`
		changequote([, ])

		if test $swig_major_ver -ge $swig_major_req &&
		   test $swig_minor_ver -ge $swig_minor_req &&
		   test $swig_micro_ver -ge $swig_micro_req
		then
			AC_MSG_RESULT(yes)
		else
			AC_MSG_RESULT(no)
			AC_MSG_ERROR(installed version of swig is too old!)
		fi
	fi
fi
])

dnl @synopsis AC_LIB_WAD
dnl
dnl This macro searches for installed WAD library.
dnl
AC_DEFUN([AC_LIB_WAD],
[
AC_ARG_ENABLE(wad,
	AC_HELP_STRING([--enable-wad], [enable wad module]),
	[
   		case "${enableval}" in
           	no)  ;;
            *)   if test "x${enableval}" = xyes;
				 then
					check_wad="yes"
				 fi ;;
        esac
	], [])

if test -n "$check_wad";
then
	# this won't work unless PYTHON_LINK and PYTHON_EXTRA_LIBS defined
	AC_CHECK_LIB(wadpy, _init, [WADPY=-lwadpy], [], $PYTHON_LINK $PYTHON_EXTRA_LIBS)
	AC_SUBST(WADPY)
fi
])

dnl ----------------------------------------------------------------------------

# # SWIG_PYTHON([use-shadow-classes])
# #
# # Checks for Python and provides the $(SWIG_PYTHON_CPPFLAGS), $(SWIG_PYTHON_LIB) and
# # $(SWIG_PYTHON_OPT) output variables.  $(SWIG_PYTHON_OPT) contains all necessary swig
# # options to generate code for Python.  Shadow classes are enabled unless the
# # value of the optional first argument is exactly 'no'.  If you need multi module
# # support use $(SWIG_PYTHON_LIB) (provided by the SWIG_MULTI_MODULE_SUPPORT() macro)
# # to link against the appropriate library.  It contains the SWIG Python runtime library
# # that is needed by the type check system for example.
# AC_DEFUN([SWIG_PYTHON],[
#         AC_REQUIRE([SWIG_PROG])
#         AC_REQUIRE([PYTHON_DEVEL])
#         if test "$SWIG" != "false" ; then
#                 AC_SUBST(SWIG_PYTHON_LIB,[-lswigpy])
#                 test ! "x$1" = "xno" && swig_shadow=" -shadow" || swig_shadow=""
#                 AC_SUBST(SWIG_PYTHON_OPT,[-python$swig_shadow])
#         fi
#         AC_SUBST(SWIG_PYTHON_CPPFLAGS,[$PYTHON_CPPFLAGS])
# ])
