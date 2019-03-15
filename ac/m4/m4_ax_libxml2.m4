AC_DEFUN([AX_LIBXML],[
    AC_SUBST(LIBXML_CFLAGS)
    AC_SUBST(LIBXML_LDFLAGS)
    AC_MSG_CHECKING([for libxml2 library])
    AC_LANG_PUSH(C++)
    save_LIBS="$LIBS"
    save_LDFLAGS="$LDFLAGS"
    LIBS="$LIBS -lxml2"
    AC_TRY_LINK([#include <libxml/xpath.h>],
                [xmlInitParser();],
                has_xml_lib=1,
                has_xml_lib=0)
    if test $has_xml_lib = 1; then
      AC_MSG_RESULT([yes])
      LIBXML_LDFLAGS="-lxml2"
    else
      AC_PATH_PROG([XML_CONFIG],[xml2-config])
      if test -z "$XML_CONFIG"; then
        AC_MSG_ERROR([Cannot find xml2-config in your system path])
      else
        cflags=`$XML_CONFIG --cflags`
        for inc in $cflags; do
            if echo $inc | egrep '/libxml2$' > /dev/null 2>&1; then
                inc=`echo $inc | sed 's|/libxml2||'`
                cflags="$cflags $inc"
            fi
        done
        libs=`$XML_CONFIG --libs`
        LIBS="$save_LIBS $libs"
        LDFLAGS="$save_LDFLAGS $cflags"
        AC_TRY_LINK([#include <libxml/xpath.h>],
                    [xmlInitParser();],
                    has_xml_lib=1,
                    has_xml_lib=0)
        if test $has_xml_lib = 1; then
          AC_MSG_RESULT([yes])
          LIBXML_LDFLAGS="$libs"
          LIBXML_CFLAGS=$cflags
        else
          AC_MSG_RESULT([no])
          AC_MSG_ERROR([Cannot find libxml2])
        fi
      fi
    fi
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS"
    AC_LANG_POP(C++)
])
