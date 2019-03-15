AC_INIT([casa], [4.6.0])
AC_COPYRIGHT([Copyright (C) 2017 Associated Universities, Inc. Washington DC, USA.])
AC_REVISION([syscmd([ac/scripts/configure.commit])])

m4_include([ac/m4/m4_ax_pkg_prog_pkg_config.m4])
m4_include([ac/m4/m4_ax_boost_base.m4])
m4_include([ac/m4/m4_ax_canonical_name.m4])
m4_include([ac/m4/m4_ax_canonical_split.m4])
m4_include([ac/m4/m4_ax_gsl.m4])
m4_include([ac/m4/m4_ax_cfitsio.m4])
m4_include([ac/m4/m4_ax_wcslib.m4])
m4_include([ac/m4/m4_ax_rpfits.m4])
m4_include([ac/m4/m4_ax_check_compile_flag.m4])
m4_include([ac/m4/m4_ax_cxx11.m4])
m4_include([ac/m4/m4_ax_cxx_compile_stdcxx_11.m4])
m4_include([ac/m4/m4_ax_expand_path.m4])
m4_include([ac/m4/m4_ax_lib_stdcxx.m4])
m4_include([ac/m4/m4_ax_libsakura.m4])
m4_include([ac/m4/m4_ax_libxml2.m4])
m4_include([ac/m4/m4_ax_openmp.m4])
m4_include([ac/m4/m4_ax_path_to_binary.m4])
m4_include([ac/m4/m4_ax_pthread.m4])
m4_include([ac/m4/m4_ax_python_devel.m4])
m4_include([ac/m4/m4_ax_python_numpy.m4])
m4_include([ac/m4/m4_ax_eigen.m4])
m4_include([ac/m4/m4_ax_swig.m4])

dnl setup cannonical host variables without requiring
dnl install.sh et al. as AC_CANONICAL_HOST does
AX_CANONICAL_NAME
AC_MSG_NOTICE([host os is ${host_osname}])

PWD=`pwd`
AC_SUBST(PWD)

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
AC_SUBST(OPENMP_CXXFLAGS)

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

AC_MSG_CHECKING([for std::defaultfloat])
AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[#include <ios>]],
                                    [[std::defaultfloat;]])],
                   [AC_MSG_RESULT([yes])],
                   [AC_MSG_RESULT([no])
                    AC_MSG_ERROR([C++ compiler is not new enough]) ] )

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
AC_MSG_NOTICE([c++ compiler is $CXX])
save_PATH_save=${PATH}
### have to protect against picking up C compiler from places other than
### where the C++ compiler is found...
temporary_script_dir=/tmp/casatools-config-$$
AC_MSG_NOTICE([ temporary directory $temporary_script_dir ])
current_directory=`pwd`
AC_MSG_NOTICE([creating $temporary_script_dir])
mkdir $temporary_script_dir
if [[ $? -ne 0 ]] ; then
   AC_MSG_NOTICE([creation of temporary directory $temporary_script_dir failed])
   exit 1
fi
cd $temporary_script_dir
ln -s `which sed`
ln -s `which rm`
ln -s `which dirname`
ln -s `which basename`
ln -s `which grep`
cd $current_directory

PATH=$(dirname ${CXX}):/bin:$temporary_script_dir

if [[ "$host_osname" = "darwin" ]]; then
    AC_PROG_CC([clang cc gcc gcc-mp-5])
else
    AC_PROG_CC([gcc cc clang])
fi
AX_PATH_TO_BINARY(${CC},CC)
AC_MSG_NOTICE([C compiler ${CC}])
PATH=$(dirname ${CXX})${PATH_SEPARATOR}${PATH}

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
_ccache_path=`perl -e "use Cwd 'abs_path'; print abs_path('${CCACHE}');"`
_cc_path=`perl -e "use Cwd 'abs_path'; print abs_path('${CC}');"`

if test "${_ccache_path}" = "${_cc_path}"; then
    ###
    ### user can utilize ccache by adding a directory of symlinks to PATH
    ###
    CCACHE=""
fi

for func in pread pwrite; do
    AC_CHECK_FUNC([$func],[],AC_MSG_ERROR([$func function is required by casacore]))
done

AX_BOOST_BASE(1.41)
AS_IF([test "x$HAVE_BOOST" = "x1"], [
    BOOST_CFLAGS="-I${BOOST_INCLUDE_PATH}"
    BOOST_LDFLAGS="-L${BOOST_LIBRARY_PATH}"
],[
    BOOST_CFLAGS=
    BOOST_LDFLAGS=
])
AC_SUBST(BOOST_CFLAGS)
AC_SUBST(BOOST_LDFLAGS)

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
AX_WCSLIB
AX_RPFITS
AX_GSL(2.0)
AX_SWIG(3.0.0)

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

AX_LIBXML
AX_EIGEN
AS_IF([ test "x$HAVE_EIGEN" = "x0" ], [
    AS_IF([ ! test -d eigen-source ], [
        AC_MSG_NOTICE([Could not find Eigen3, fetching source code...])
        ( git clone https://github.com/eigenteam/eigen-git-mirror.git eigen-source || rm -rf eigen-source ) > /dev/null 2>&1
        pushd eigen-source > /dev/null 2>&1
        git checkout > /dev/null 2>&1
        git checkout branches/3.2 > /dev/null 2>&1
        popd > /dev/null 2>&1
    ],[
        AC_MSG_NOTICE([Could not find Eigen3 but source code already fetched...])
    ])
    HAVE_EIGEN=1
    EIGEN_CFLAGS="-I`pwd`/eigen-source"
    FETCHED_EIGEN=1
],[
    FETCHED_EIGEN=0
])
AC_SUBST(FETCHED_EIGEN)

if [[ -d "$temporary_script_dir" ]] ; then
    AC_MSG_NOTICE([removing $temporary_script_dir])
    rm -rf $temporary_script_dir
fi

AC_CONFIG_FILES(build.properties:ac/templates/build.in)
AC_CONFIG_FILES(binding/source/version.cc:casa-source/code/stdcasa/version.cc.in)
AC_CONFIG_FILES([setup.py:ac/templates/setup.py.in],[chmod +x setup.py])
AC_OUTPUT( )
