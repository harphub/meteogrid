AC_INIT(meteogrid, 3.7.2.9000, alex.deckmyn@meteo.be)

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
  LIBS="-L${proj_lib_path} ${LIBS}"
else
  if test [ -n "${PROJ_LIBS}" ] ; then
    LIBS="-L${PROJ_LIBS} ${LIBS}"
  fi
fi

proj4ok=yes
AC_CHECK_HEADERS(proj_api.h,,proj4ok=no)
# don't fail when proj is not found, maybe it's OK for R (Makevars)
#if test "${proj4ok}" = no; then
#   AC_MSG_ERROR([proj_api.h not found.
#*** Install PROJ.4 and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly.])
#fi
AC_CHECK_LIB(proj,pj_init,,proj4ok=no)
#if test "${proj4ok}" = no; then
#   AC_MSG_ERROR([libproj not found.
#*** Install PROJ.4 and if necessary set PKG_CPPFLAGS/PKG_LIBS accordingly.])
#fi

# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
  exit 1
fi
RBIN="${R_HOME}/bin/R"

# pick all flags for testing from R
CC=`"${RBIN}" CMD config CC`
CPP=`"${RBIN}" CMD config CPP`
CFLAGS=`"${RBIN}" CMD config CFLAGS`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
LDFLAGS=`"${RBIN}" CMD config LDFLAGS`
AC_PROG_CC
AC_PROG_CPP

CPPFLAGS="${CPPFLAGS} ${PROJ_CPPFLAGS}"

AC_SUBST(PROJ_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT



