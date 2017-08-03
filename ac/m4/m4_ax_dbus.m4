AC_DEFUN([AX_DBUS],[
    AC_SUBST(DBUS_CFLAGS)
    AC_SUBST(DBUS_LDFLAGS)
    AC_MSG_CHECKING([for dbus library])
    AC_LANG_PUSH(C++)
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -ldbus-1"
    AC_TRY_LINK([#include <dbus/dbus.h>], 
                [DBusError err; dbus_error_init(&err);],
                has_dbus_lib=1,
                has_dbus_lib=0)
    if test $has_dbus_lib = 1; then
      AC_MSG_RESULT([yes])
      DBUS_LDFLAGS="-ldbus-1"
    else
      if test -d /usr/include/dbus-1.0 -a -d /usr/lib64/dbus-1.0/include; then
        LDFLAGS="$LDFLAGS -I/usr/include/dbus-1.0 -I/usr/lib64/dbus-1.0/include"
        AC_TRY_LINK([#include <dbus/dbus.h>], 
                    [DBusError err; dbus_error_init(&err);],
                    has_dbus_lib=1,
                    has_dbus_lib=0)
        if test $has_dbus_lib = 1; then
          AC_MSG_RESULT([yes])
          DBUS_CFLAGS="-I/usr/include/dbus-1.0 -I/usr/lib64/dbus-1.0/include"
          DBUS_LDFLAGS="-ldbus-1"
        else
          AC_MSG_RESULT([no])
          AC_MSG_ERROR([Cannot find dbus])
        fi
      fi
      LDFLAGS="$save_LDFLAGS"
      LIBS="$save_LIBS"
      AC_LANG_POP(C++)
    fi
])
