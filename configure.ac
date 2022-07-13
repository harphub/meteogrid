AC_INIT(meteogrid, 3.9.2.9001, alex.deckmyn@meteo.be)

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
# FIRST: both combined in a single --with-proj statement:
AC_ARG_WITH([proj],
            AC_HELP_STRING([--with-proj=PROJ_PATH],
                           [the location of proj]),
            [proj_path=$withval])
if test [ -n "$proj_path" ] ; then
  AC_MSG_NOTICE([Using proj_path=${proj_path}])
  LIBS="-L${proj_path}/lib -Wl,-rpath,${proj_path}/lib ${LIBS}"
  PROJ_CPPFLAGS="-I. -I${proj_path}/include"
elif test [ -n "${PROJ_PATH}" ] ; then
    AC_MSG_NOTICE([Using PROJ_PATH=${PROJ_PATH}])
    LIBS="-L${PROJ_PATH}/lib -Wl,-rpath,${PROJ_PATH}/lib ${LIBS}"
    PROJ_CPPFLAGS="-I. -I${PROJ_PATH}/include"
else
  # SECOND: separate --with-proj-[lib|include] statements:
  AC_ARG_WITH([proj-include],
            AC_HELP_STRING([--with-proj-include=INCLUDE_PATH],
                           [the location of the proj header files]),
            [proj_include_path=$withval])
  PROJ_CPPFLAGS="-I."
  if test [ -n "$proj_include_path" ] ; then
    AC_MSG_NOTICE([Using proj_include_path=${proj_include_path}])
    PROJ_CPPFLAGS="-I. -I${proj_include_path}"
  elif test [ -n "$PROJ_INCLUDE_PATH" ] ; then
    AC_MSG_NOTICE([Using PROJ_INCLUDE_PATH=${PROJ_INCLUDE_PATH}])
    PROJ_CPPFLAGS="-I. -I${PROJ_INCLUDE_PATH}"
  fi

  AC_ARG_WITH([proj-lib],
            AC_HELP_STRING([--with-proj-lib=LIB_PATH],
                           [the location of the proj libraries]),
            [proj_lib_path=$withval])
  if test [ -n "$proj_lib_path" ] ; then
    AC_MSG_NOTICE([Using proj_lib_path=${proj_lib_path}])
    LIBS="-L${proj_lib_path} -Wl,-rpath,${proj_lib_path} ${LIBS}"
  elif test [ -n "${PROJ_LIB_PATH}" ] ; then
    AC_MSG_NOTICE([Using PROJ_LIB_PATH=${PROJ_LIB_PATH}])
    LIBS="${PROJ_LIB_PATH} -Wl,-rpath,${PROJ_LIB_PATH} ${LIBS}"
  fi
fi
CPPFLAGS="${PROJ_CPPFLAGS} ${CPPFLAGS}"

my_error="*** Install PROJ or add --with-proj=<path/to/proj>
  *** or --with-proj-include=<.../include> --with-proj-lib=<.../lib>
  *** or set variable proj_path (proj_include_path and proj_lib_path)
  *** or (cfr ECMWF) set variable PROJ_PATH (PROJ_INCLUDE_PATH and PROJ_LIB_PATH)
"

AC_CHECK_HEADERS(proj.h proj_api.h,  , )
# first check for the new proj.h interface (>=v5.0)
if test ${ac_cv_header_proj_h} = yes ; then
  AC_MSG_NOTICE([Using proj.h])
  AC_CHECK_LIB(proj, proj_create, , AC_MSG_ERROR([Can not compile with proj.h
  *** Make sure libproj is found.
  ${my_error}]))
elif test ${ac_cv_header_proj_api_h} = yes ; then
  # if we use the old interface, we need to add -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H
  # also, AC_CHECK_LIB must then run with this setting
  AC_MSG_NOTICE([Using proj_api.h (deprecated!)])
  PROJ_CPPFLAGS="${PROJ_CPPFLAGS} -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H"
  AC_DEFINE(ACCEPT_USE_OF_DEPRECATED_PROJ_API_H)
  AC_CHECK_LIB(proj, pj_init, , AC_MSG_ERROR([Can not compile with proj_api.h
  *** Make sure libproj is found.
  ${my_error}]))
else
  AC_MSG_ERROR([proj.h and proj_api.h both not found.
  ${my_error}])
fi

AC_SUBST(PROJ_CPPFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT

