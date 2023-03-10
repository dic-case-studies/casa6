AC_INIT([casa], [4.6.0])
AC_COPYRIGHT([Copyright (C) 2017-2021 Associated Universities, Inc. Washington DC, USA.])
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
m4_include([ac/m4/ax_prog_cxx_mpi.m4])
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

WITH_GRPC=1
AC_ARG_ENABLE( [system-grpc], AS_HELP_STRING([--enable-system-grpc],[use system grpc libraries]),
               [SYSTEM_GRPC=1], [SYSTEM_GRPC=0] )

if [[ "x${CXX}" != "x" ]]; then
    dnl#########################################################################
    dnl##### set CXX_OVERRIDE to prevent all of our C++ machinations from  #####
    dnl##### selecting some other C++ compiler...                          #####
    dnl#####                                                               #####
    dnl##### since compilers have different ABIs and different STLs        #####
    dnl#########################################################################
    CXX_OVERRIDE="${CXX}"
else
    CXX_OVERRIDE=""
fi

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

if [[ ! -e "${CXX}" ]]; then
    CXX=`which ${CXX}`
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


dnl Add --with-mpi option to configure
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. Auto: if no MPI compiler/libraries are
    found, MPI is not used. Default: auto])
    ],,
    [with_mpi=auto])

dnl Try to find an MPI compiler wrapper. Far from ideal.
dnl AC_PROG_CXX_MPI is meant to be used intead of AC_PROG_CXX for MPI applications.
dnl But we have "m4_ax_cxx11" fiddling with CXX etc variables up here.
if [test x"$with_mpi" != xno]; then
   _PRE_MPI_CXX=$CXX
   _PRE_MPI_CC=$CC
   _PRE_MPI_PATH=$PATH
   _CXX_PATH=$(dirname ${CXX})
   _PATH_CASA_OPENMPI=$PATH${PATH_SEPARATOR}/opt/casa/03/bin/
   AC_MSG_NOTICE([Looking for MPI compiler wrappers in extended path: ${_PATH_CASA_OPENMPI}])
   dnl want to find this or alternatives: CXX=/opt/casa/03/bin/mpic++
   AC_PATH_PROGS(CASA_MPI_CXX, [mpiCC, mpicxx, mpic++],,${_PATH_CASA_OPENMPI})
   AC_MSG_NOTICE([After looking for an MPI wrapper, the C++ compiler is: ${CASA_MPI_CXX}])
   CXX=$CASA_MPI_CXX
   PATH=${_CXX_PATH}${PATH_SEPARATOR}${PATH}
   dnl want to find this or alternatives: CC=/opt/casa/03/bin/mpicc
   AC_PATH_PROGS(CASA_MPI_CC, [mpicc],,${_PATH_CASA_OPENMPI})
   AC_MSG_NOTICE([After looking for an MPI wrapper, the C compiler is: ${CASA_MPI_CC}])
   CC=$CASA_MPI_CC
   AX_PROG_CXX_MPI(true,
     [
     CASA_HAVE_MPI=1
     AC_SUBST(CASA_HAVE_MPI)
     ],[
     CXX=$_PRE_MPI_CXX
     CC=$_PRE_MPI_CC
     PATH=$_PRE_MPI_PATH
     if test x"$with_mpi" = xyes; then
       AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
     else
       AC_MSG_WARN([No working MPI compiler could be found, won't use MPI.])
     fi
   ])
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
ln -s `which sort`
ln -s `which printf`
ln -s `which cat`
cd $current_directory

PATH=$(dirname ${CXX}):/bin:$temporary_script_dir

if [[ "$host_osname" = "darwin" ]]; then
    AC_PROG_CC([clang cc gcc gcc-mp-5])
else
    AC_PROG_CC([gcc cc clang])
fi
AX_PATH_TO_BINARY(${CC},CC)
AC_MSG_NOTICE([C compiler ${CC}])

AC_LANG_PUSH([C])
AX_CHECK_COMPILE_FLAG(-fcx-fortran-rules,[CC_FORTRAN_COMPLEX=-fcx-fortran-rules],[CC_FORTRAN_COMPLEX=])
AC_SUBST(CC_FORTRAN_COMPLEX)
AC_LANG_POP([C])

for pp in $(dirname ${CXX})${PATH_SEPARATOR}${PATH} ${save_PATH_save}; do
    PATH=${pp}:${temporary_script_dir}
    AC_PROG_FC([gfortran-mp-9 gfortran-mp-8 gfortran-mp-5 gfortran-mp-4.9 gfortran-mp-4.8 gfortran])
    AX_PATH_TO_BINARY(${FC},FC)
    if [[ -f "$FC" -a -x "$FC" ]]; then
        AC_MSG_NOTICE([Fortran compiler ${FC}])
        break
    else
        unset ac_cv_prog_FC
        unset FC
    fi
done

if [[ ! -f "$FC" -o ! -x "$FC" ]]; then
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

if [[ "${SYSTEM_GRPC}" == "1" ]]; then
    PKG_CHECK_MODULES( [GRPC], [grpc >= 6.0.0 grpc++ >= 1.16.1 protobuf >= 3.6.1 ], [
                       AC_CHECK_PROGS([PROTOC], [protoc], [/usr/bin/false])
                       if test "$PROTOC" = "/usr/bin/false"; then
                           AC_MSG_FAILURE([grpc requested but protoc not found.])
                       else
                           PROTOC=`which protoc`
                       fi
                       AC_CHECK_PROGS([GRPC_CPP_PLUGIN], [grpc_cpp_plugin], [/usr/bin/false])
                       if test "$GRPC_CPP_PLUGIN" = "/usr/bin/false"; then
                           AC_MSG_FAILURE([grpc requested but grpc plugin for C++ not found.])
		       else
		           GRPC_CPP_PLUGIN=`which grpc_cpp_plugin`
                       fi ] )
    GRPC_DO_BUILD=0
else
    GRPC_DO_BUILD=1
    PROTOC="/usr/bin/false"
    GRPC_CPP_PLUGIN="/usr/bin/false"
fi
AC_SUBST(GRPC_CFLAGS)
AC_SUBST(GRPC_LIBS)
AC_SUBST(GRPC_DO_BUILD)
AC_SUBST(PROTOC)

#AC_MSG_CHECKING([for casa version number])
#read -a full_version <<< "`scripts/version | perl -pe 's|^(\d+).(\d+).*?$|$1 $2 0 0|'`"
#CASA_VERSION_MAJOR=${full_version[[0]]}
#CASA_VERSION_MINOR=${full_version[[1]]}
#CASA_VERSION_PATCH=${full_version[[2]]}
#CASA_VERSION_FEATURE=${full_version[[3]]}
#CASA_VERSION_DESC="CASAtools:v1.0.0"
#CASA_VERSION="${CASA_VERSION_MAJOR}.${CASA_VERSION_MINOR}.${CASA_VERSION_PATCH}-${CASA_VERSION_FEATURE}"
#CASA_VARIANT="none"
#AC_MSG_RESULT($CASA_VERSION)


#AC_SUBST(CASA_VERSION_MAJOR)
#AC_SUBST(CASA_VERSION_MINOR)
#AC_SUBST(CASA_VERSION_PATCH)
#AC_SUBST(CASA_VERSION_FEATURE)
#AC_SUBST(CASA_VERSION_DESC)
#AC_SUBST(CASA_VARIANT)
#AC_SUBST(CASA_VERSION)

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
#AC_CONFIG_FILES(binding/source/version.cc:casa-source/code/stdcasa/version.cc.in)
AC_CONFIG_FILES([setup.py:ac/templates/setup.py.in],[chmod +x setup.py])
AC_OUTPUT( )
