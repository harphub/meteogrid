AC_INIT(meteogrid, 3.7.2.9000, alex.deckmyn@meteo.be)

# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
  exit 1
fi

# pick all flags for testing from R (e.g. use ~/.R/Makevars)
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
LDFLAGS=`"${RBIN}" CMD config LDFLAGS`
AC_PROG_CC
AC_PROG_CPP

# check for user-specified PROJ
AC_ARG_WITH([proj-include],
            AC_HELP_STRING([--with-proj-include=INCLUDE_PATH],
                           [the location of the proj header files]),
            [proj_include_path=$withval])
PROJ_CPPFLAGS="-I."
if test [ -n "$proj_include_path" ] ; then
  PROJ_CPPFLAGS="-I. -I${proj_include_path}"
else
  if test [ -n "$PROJ_INCLUDE" ] ; then
    PROJ_CPPFLAGS="-I. -I${PROJ_INCLUDE}"
  fi
fi

AC_ARG_WITH([proj-lib],
            AC_HELP_STRING([--with-proj-lib=LIB_PATH],
                           [the location of the proj libraries]),
            [proj_lib_path=$withval])
if test [ -n "$proj_lib_path" ] ; then
  LIBS="-L${proj_lib_path} -Wl,-rpath,${proj_lib_path} ${LIBS}"
else
  if test [ -n "${PROJ_LIBS}" ] ; then
    LIBS="-L${PROJ_LIBS} -Wl,-rpath,${PROJ_LIBS} ${LIBS}"
  fi
fi

# both combined in a single --with-proj statement:
AC_ARG_WITH([proj],
            AC_HELP_STRING([--with-proj=PROJ_PATH],
                           [the location of proj]),
            [proj_path=$withval])
if test [ -n "$proj_path" ] ; then
  LIBS="-L${proj_path}/lib -Wl,-rpath,${proj_path}/lib ${LIBS}"
  PROJ_CPPFLAGS="-I. -I${proj_path}/include"
else
  if test [ -n "${PROJ_DIR}" ] ; then
    LIBS="-L${PROJ_DIR}/lib -Wl,-rpath,${PROJ_DIR}/lib ${LIBS}"
    PROJ_CPPFLAGS="-I. -I${PROJ_DIR}/include"
  fi
fi

# first check for the new proj.h interface (>=v5.0)
AC_CHECK_HEADERS(proj.h, proj_ok=yes, proj_ok=no)
if test "${proj_ok}" = yes; then
#  PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DPROJ5"
  AC_CHECK_LIB(proj, proj_create, , proj_ok=no)
else
# if we use the old interface, we need to add -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H
# also, AC_CHECK_LIB must then run with this setting
  PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H"
  AC_DEFINE(ACCEPT_USE_OF_DEPRECATED_PROJ_API_H)
  AC_CHECK_HEADERS(proj_api.h, proj_ok=yes, proj_ok=no)
  AC_CHECK_LIB(proj, pj_init, , proj_ok=no)
fi

if test "${proj_ok}" = no; then
   AC_MSG_ERROR([proj.h and proj_api.h both not found.
*** Install PROJ or add --with-proj=<path/to/proj>
*** or --with-proj-include=<path/to/proj/include> --with-proj-lib=<path/to/proj/lib>
*** or set environment variables PROJ_DIR or PROJ_INCLUDE and PROJ_LIBS])
fi

CPPFLAGS="${CPPFLAGS} ${PROJ_CPPFLAGS}"

AC_SUBST(PROJ_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT

