# AX_CANONICAL_SPLIT(THING,CANNONICAL_STRING)
# --------------------------
# Generate the variables THING, THING_{alias cpu vendor os}.
AC_DEFUN([AX_CANONICAL_SPLIT],
[case $2 in
*-*-*) ;;
*) AC_MSG_ERROR([invalid value of canonical $2]);;
esac
AC_SUBST([$1], [$2])dnl
ac_save_IFS=$IFS; IFS='-'
set x $2
shift
AC_SUBST([$1_cpu], [$[1]])dnl
AC_SUBST([$1_vendor], [$[2]])dnl

shift; shift
[# Remember, the first character of IFS is used to create $]*,
# except with old shells:
$1_os=$[*]
IFS=$ac_save_IFS
case $$1_os in *\ *) $1_os=`echo "$$1_os" | sed 's/ /-/g'`;; esac
AC_SUBST([$1_os])dnl
]

[case $$1_os in
linux*) $1_osname=linux;;
darwin*) $1_osname=darwin;;
*) $1_osname=$$1_os;;
esac
AC_SUBST([$1_osname])
])
