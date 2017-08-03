AC_DEFUN([AX_LIBSAKURA],[
    AC_SUBST(LIBSAKURA_CFLAGS)
    AC_SUBST(LIBSAKURA_LDFLAGS)
    AC_MSG_CHECKING([for libsakura])
    AC_LANG_PUSH(C++)
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lsakura"
    AC_TRY_LINK([#include <libsakura/sakura.h>], 
                [LIBSAKURA_SYMBOL(Initialize)(0,0);],
                has_libsakura=1,
                has_libsakura=0)
    if test $has_libsakura = 1; then
      AC_MSG_RESULT([yes])
      LIBSAKURA_LDFLAGS="-lsakura"
    else
      for path in '/opt/casa/02/lib/libsakura/default' ; do
        LDFLAGS="$LDFLAGS -I$path/include -L$path/lib"
        AC_TRY_LINK([#include <libsakura/sakura.h>], 
                    [LIBSAKURA_SYMBOL(Initialize)(0,0);],
                    has_libsakura=1,
                    has_libsakura=0)
        if test $has_libsakura = 1; then
          AC_MSG_RESULT([yes])
          LIBSAKURA_CFLAGS="-I$path/include"
          LIBSAKURA_LDFLAGS="-L$path/lib -lsakura"
        else
          AC_MSG_RESULT([no])
        fi
      done
    fi
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS"
    AC_LANG_POP(C++)
])
