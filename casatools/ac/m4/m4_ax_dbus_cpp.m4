AC_DEFUN([AX_DBUSCPP],[
	AC_PATH_PROG([DBUSCPP_XML2CPP],[dbuspp-xml2cpp])
	if test -z "$DBUSCPP_XML2CPP"; then
	   AC_MSG_ERROR([Cannot find dbuspp-xml2cpp in your system path])
	fi
    AC_SUBST(DBUSCPP_XML2CPP)
    AC_SUBST(DBUSCPP_CFLAGS)
    AC_SUBST(DBUSCPP_LDFLAGS)
    AC_SUBST(DBUSCPP_NAME)
    AC_MSG_CHECKING([for casa-dbus-cpp library])
    AC_LANG_PUSH(C++)
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lcasa-dbus-cpp"
    AC_TRY_LINK([#include <dbus-cpp/dbus.h>],
                [DBus::Connection::SessionBus( );],
                has_dbus_cpp_lib=1,
                has_dbus_cpp_lib=0)
    if test $has_dbus_cpp_lib = 1; then
      AC_MSG_RESULT([yes])
      DBUSCPP_LDFLAGS="-lcasa-dbus-cpp"
      DBUSCPP_NAME="casa-dbus-cpp"
    else
      dbus_root=`cd $(dirname $DBUSCPP_XML2CPP)/.. && pwd`
      LDFLAGS="$LDFLAGS -I$dbus_root/include -L$dbus_root/lib"
      AC_TRY_LINK([#include <dbus-cpp/dbus.h>],
                  [DBus::Connection::SessionBus( );],
                  has_dbus_cpp_lib=1,
                  has_dbus_cpp_lib=0)
      if test $has_dbus_cpp_lib = 1; then
        AC_MSG_RESULT([yes])
        DBUSCPP_CFLAGS="-I$dbus_root/include"
        DBUSCPP_LDFLAGS="-L$dbus_root/lib -lcasa-dbus-cpp"
        DBUSCPP_NAME="casa-dbus-cpp"
      else
        LIBS="$save_LIBS -ldbus-cpp"
        AC_TRY_LINK([#include <dbus-cpp/dbus.h>],
                    [DBus::Connection::SessionBus( );],
                    has_dbus_cpp_lib=1,
                    has_dbus_cpp_lib=0)
        if test $has_dbus_cpp_lib = 1; then
          AC_MSG_RESULT([yes])
          DBUSCPP_CFLAGS="-I$dbus_root/include"
          DBUSCPP_LDFLAGS="-L$dbus_root/lib -ldbus-cpp"
          DBUSCPP_NAME="dbus-cpp"
        else
          AC_MSG_RESULT([no])
        fi
      fi
    fi
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS"
    AC_LANG_POP(C++)
])
