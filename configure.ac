AC_INIT([casa], [4.6.0])
AC_COPYRIGHT([Copyright (C) 2017 Associated Universities, Inc. Washington DC, USA.])
AC_REVISION([syscmd([ac/scripts/configure.commit])])

m4_include([ac/m4/m4_ax_pkg_prog_pkg_config.m4])
m4_include([ac/m4/m4_ax_boost_base.m4])
m4_include([ac/m4/m4_ax_canonical_name.m4])
m4_include([ac/m4/m4_ax_canonical_split.m4])
m4_include([ac/m4/m4_ax_cfitsio.m4])
m4_include([ac/m4/m4_ax_check_compile_flag.m4])
m4_include([ac/m4/m4_ax_cxx11.m4])
m4_include([ac/m4/m4_ax_cxx_compile_stdcxx_11.m4])
m4_include([ac/m4/m4_ax_dbus.m4])
m4_include([ac/m4/m4_ax_dbus_cpp.m4])
m4_include([ac/m4/m4_ax_expand_path.m4])
m4_include([ac/m4/m4_ax_lib_stdcxx.m4])
m4_include([ac/m4/m4_ax_libsakura.m4])
m4_include([ac/m4/m4_ax_libxml2.m4])
m4_include([ac/m4/m4_ax_path_to_binary.m4])
m4_include([ac/m4/m4_ax_pthread.m4])
m4_include([ac/m4/m4_ax_python_devel.m4])
m4_include([ac/m4/m4_ax_python_numpy.m4])

dnl setup cannonical host variables without requiring
dnl install.sh et al. as AC_CANONICAL_HOST does
AX_CANONICAL_NAME
AC_MSG_NOTICE([host os is ${host_osname}])

PWD=`pwd`
AC_SUBST(PWD)

AC_ARG_ENABLE([dbus], [  --enable-dbus           compile in support for DBus (deprecated feature)])
AC_ARG_ENABLE([grpc], [  --disable-grpc          compile in support for gRPC],
  [case "${enableval}" in
     yes ) WITH_GRPC=1 ;;
     no ) WITH_GRPC=0 ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-grpc) ;;
   esac],
  [WITH_GRPC=1]
)
AC_SUBST(WITH_GRPC)


AC_LANG(C++)

AX_CXX11($PATH)

if [[ "${ac_success}" = "no" ]]; then
    dnl#############################################################################################
    dnl###  could not find a C++ compiler that supported C++11 in user's PATH                    ###
    dnl#############################################################################################
    case $host_osname in
    linux)
        AX_CXX11_CLEAR_CACHE
        AC_MSG_NOTICE([could not find a C++ compiler that supports C++11, in your PATH... trying other possible linux paths...])
        AX_CXX11([/opt/rh/devtoolset-3/root/usr/bin:/opt/rh/devtoolset-2/root/usr/bin:/usr/local/bin:/opt/local/bin:/usr/bin:/bin])
        ;;
    darwin)
        AX_CXX11_CLEAR_CACHE
        AC_MSG_NOTICE([could not find a C++ compiler that supports C++11, in your PATH... trying other possible osx paths...])
        AX_CXX11([/usr/bin:/opt/local/bin:/opt/macports/bin:/bin])
        ;;
    *)
        AC_MSG_ERROR([cannot find a C++ compiler that supports C++11, please include the path to one in your PATH environment variable])
        ;;
    esac

    if [[ "${ac_success}" = "no" ]]; then
        AC_MSG_ERROR([cannot find a C++ compiler that supports C++11, please include the path to one in your PATH environment variable])
    fi
fi

if [[ -x "${CXX}" ]]; then
    AC_MSG_NOTICE([C++11 compiler ${CXX}])
else
    AC_MSG_ERROR([failed to find viable C++11 compiler, please include the path to one in your PATH environment variable])
fi

dnl# use faster fortran rules for complex operations, removes restoring complex infinities
dnl# if naive computation results in NAN + NAN * I  Handling complex multiplication and
dnl# division with correct treating of complex infinities (one element Inf regardless of
dnl# the other) according to the C is complicated, e.g
dnl#
dnl#     if a = NaN + 1e30 i; a * a is not NaN but a complex infinity (-Inf - NaN).
dnl#
dnl# Treating this situation correctly has large performance impact. In GCC's implementation
dnl# it is about 4 times slower than the naive implementation, with vectorization enabled the
dnl# impact is even larger. As correct treatment of complex infinities when NaN appear in
dnl# results is seldom accounted for, or not required and most other languages do not have
dnl# these rules, the correct treatmeant can be disabled with the -fcx-fortran-rules flag.
dnl# This changes the semantics to those of the FORTRAN language which is removes the need for
dnl# rescuing the result when NaN appear. Python also follows FORTRAN rules.
dnl#
dnl# Additionally the correct behavior is not implemented in all compilers, most notably clang
dnl# which is the default compiler on MacOS. So turning off correct treatment with GCC does
dnl# not only make our code faster but also behave the same on more compilers.
dnl#
AX_CHECK_COMPILE_FLAG(-fcx-fortran-rules,[CXX_FORTRAN_COMPLEX=-fcx-fortran-rules],[CXX_FORTRAN_COMPLEX=])
AC_SUBST(CXX_FORTRAN_COMPLEX)
AX_LIB_STDCXX

dnl#################################################################################################
dnl###  add C++ compiler path to start of PATH to try and find matching c/fortran compilers      ###
dnl#################################################################################################
save_PATH_save=${PATH}
PATH=$(dirname $CXX)${PATH_SEPARATOR}${PATH}

AC_PROG_CC
AX_PATH_TO_BINARY(${CC},CC)
AC_MSG_NOTICE([C compiler ${CC}])

AC_LANG_PUSH([C])
AX_CHECK_COMPILE_FLAG(-fcx-fortran-rules,[CC_FORTRAN_COMPLEX=-fcx-fortran-rules],[CC_FORTRAN_COMPLEX=])
AC_SUBST(CC_FORTRAN_COMPLEX)
AC_LANG_POP([C])

AC_PROG_FC([gfortran-mp-5 gfortran-mp-4.9 gfortran-mp-4.8 gfortran])
AX_PATH_TO_BINARY(${FC},FC)
AC_MSG_NOTICE([Fortran compiler ${FC}])
if test -f "$FC" -a ! -x "$FC"; then
   echo "could not find fortran compiler..."
   exit 1
fi
PATH=${save_PATH_save}

AC_PATH_PROG([CCACHE],[ccache])

for func in pread pwrite; do
    AC_CHECK_FUNC([$func],[],AC_MSG_ERROR([$func function is required by casacore]))
done

AX_BOOST_BASE(1.41)
AC_SUBST(BOOST_INCLUDE_PATH)

AX_PYTHON_DEVEL(>= '2.7.0')
dnl#AX_BOOST_PYTHON
AX_PYTHON_NUMPY
AC_CHECK_LIB(readline,readline)

AC_PROG_LEX
if [[ "x$LEX" != xflex ]]; then
    AC_MSG_ERROR([cannot find a version of flex to use, please adjust your PATH])
fi
AX_PATH_TO_BINARY(${LEX},FLEX)

AC_PROG_YACC
if [[ "$YACC" != "bison -y" ]]; then
    AC_MSG_ERROR([cannot find a version of bison to use, please adjust your PATH])
fi
AX_PATH_TO_BINARY(bison,BISON)

AX_PTHREAD
AX_CFITSIO

AC_MSG_CHECKING([for casa version number])
read -a full_version <<< "`(cd casa-source && code/install/resolvegitrevision.sh --pretty-print) | perl -pe 's|^(\d+)\.(\d+)\.(\d+).*?(\d+).?$|$1 $2 $3 $4|'`"
CASA_VERSION_MAJOR=${full_version[[0]]}
CASA_VERSION_MINOR=${full_version[[1]]}
CASA_VERSION_PATCH=${full_version[[2]]}
CASA_VERSION_FEATURE=${full_version[[3]]}
CASA_VERSION_DESC="CASAtools:v1.0.0"
CASA_VERSION="${CASA_VERSION_MAJOR}.${CASA_VERSION_MINOR}.${CASA_VERSION_PATCH}-${CASA_VERSION_FEATURE}"
AC_MSG_RESULT($CASA_VERSION)


AC_SUBST(CASA_VERSION_MAJOR)
AC_SUBST(CASA_VERSION_MINOR)
AC_SUBST(CASA_VERSION_PATCH)
AC_SUBST(CASA_VERSION_FEATURE)
AC_SUBST(CASA_VERSION_DESC)
AC_SUBST(CASA_VERSION)

AS_IF([test "x$enable_dbus" = "xyes"], [
    HAVE_DBUS=1
    AX_DBUS
    AX_DBUSCPP
],[
    HAVE_DBUS=0
    AC_SUBST(HAVE_DBUS)
    AC_SUBST(DBUSCPP_XML2CPP)
    AC_SUBST(DBUSCPP_CFLAGS)
    AC_SUBST(DBUSCPP_NAME)
    AC_SUBST(DBUSCPP_LDFLAGS)
    DBUS_CFLAGS="-DWITHOUT_DBUS"
    AC_SUBST(DBUS_CFLAGS)
    AC_SUBST(DBUS_LDFLAGS)
])
AC_SUBST(HAVE_DBUS)

AX_LIBXML
AX_LIBSAKURA


AC_CONFIG_FILES(build.properties:ac/templates/build.in)
AC_CONFIG_FILES(binding/source/version.cc:casa-source/code/stdcasa/version.cc.in)
AC_CONFIG_FILES([setup.py:ac/templates/setup.py.in],[chmod +x setup.py])
AC_OUTPUT( )
