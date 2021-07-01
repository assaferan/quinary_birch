AC_INIT([birch],[2.0.0-beta.3])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign])
AC_PROG_CXX
AM_PROG_AR
AC_PROG_LIBTOOL
LT_INIT

# AC_CANONICAL_HOST is needed to access the 'host_os' variable    
AC_CANONICAL_HOST

build_linux=no
build_windows=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=yes
        ;;
    cygwin*|mingw*)
        build_windows=yes
        ;;
    darwin*)
        build_mac=yes
        ;;
    *)
        AC_MSG_ERROR(["OS $host_os is not supported"])
        ;;
esac

AM_CONDITIONAL([LINUX], [test "$build_linux" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])
AM_CONDITIONAL([OSX], [test "$build_mac" = "yes"])

AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT(Makefile src/Makefile)
AM_PROG_CC_C_O
AC_SEARCH_LIBS([m], [gmp], [gmpxx])
