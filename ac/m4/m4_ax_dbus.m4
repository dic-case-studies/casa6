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
      AC_MSG_RESULT([continuing])
      for root in [/opt/casa/03 /opt/casa/02 /usr /opt/local]; do
        for lib in [lib64 lib]; do
          if test -d "$root/include/dbus-1.0" -a -d "$root/$lib/dbus-1.0/include"; then
            AC_MSG_CHECKING([for dbus in $root])
            LDFLAGS="$save_LDFLAGS -L$root/$lib -I$root/include/dbus-1.0 -I$root/$lib/dbus-1.0/include"
            AC_TRY_LINK([#include <dbus/dbus.h>],
                        [DBusError err; dbus_error_init(&err);],
                        has_dbus_lib=1,
                        has_dbus_lib=0)
            if test $has_dbus_lib = 1; then
              AC_MSG_RESULT([yes])
              DBUS_CFLAGS="-I$root/include/dbus-1.0 -I$root/$lib/dbus-1.0/include"
              DBUS_LDFLAGS="-L$root/$lib -ldbus-1"
              found_dbus="yes"
              break
            else
              AC_MSG_RESULT([no])
            fi
          fi
        done
        if test $has_dbus_lib = 1; then
          break
        fi
      done

      if test $has_dbus_lib != 1; then
        AC_MSG_ERROR([Cannot find dbus])
      fi
      LDFLAGS="$save_LDFLAGS"
      LIBS="$save_LIBS"
      AC_LANG_POP(C++)
    fi
])
