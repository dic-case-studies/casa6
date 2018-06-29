AC_DEFUN([AX_EIGEN], [

  PKG_CHECK_MODULES(EIGEN, eigen3, [ax_eigen_ok=yes], [ax_eigen_ok=no])

  if test x"$ax_eigen_ok" = x"no"; then
    AC_LANG_SAVE
    AC_LANG([C++])

    EIGEN_CFLAGS=
    AC_CHECK_HEADER(Eigen/Core,[
      AC_MSG_CHECKING([if Eigen library is usable])
      AC_LINK_IFELSE(
       [AC_LANG_PROGRAM(
        [[#include <Eigen/Dense>
          #include <Eigen/Eigenvalues>
        ]],
        [[
          typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
          Matrix h; h << 2.0, 1.0, 1.0, 2.0;
          Matrix s; s << 1.0, 0.0, 0.0, 1.0;
          Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(h, s);
          auto evals = gen_eig_solver.eigenvalues();
         ]]
        )
       ],[
          ax_eigen_ok=yes
          EIGEN_CFLAGS=
         ]
      )
      AC_MSG_RESULT([$ax_eigen_ok])
    ])
    AC_LANG_RESTORE
  fi
  AS_IF([test "x$ax_eigen_ok" = "xyes"], [
      HAVE_EIGEN=1
  ],[
      HAVE_EIGEN=0
  ])
  AC_SUBST(HAVE_EIGEN)
  AC_SUBST(EIGEN_CFLAGS)
])
