CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)

set(AVERAGE
average.f90
)

set(RUNOFF
move_runoff.f90
)

if (CMAKE_Fortran_COMPILER MATCHES "pgf95$")
        set (F2PY_COMPILER "pg")
        set (F2PY_F90_FLAGS "-r8")
elseif (CMAKE_Fortran_COMPILER MATCHES "pgf90$")
        set (F2PY_COMPILER "pg")
        set (F2PY_F90_FLAGS "-r8")
elseif (CMAKE_Fortran_COMPILER MATCHES "gfortran$")
        set (F2PY_COMPILER "gnu95")
        set (F2PY_F90_FLAGS "-fdefault-real-8")
endif(CMAKE_Fortran_COMPILER MATCHES "pgf95$")

add_custom_target(pyroms_toolbox
COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${F2PY_F90_FLAGS} -m _average -c ${AVERAGE}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${F2PY_F90_FLAGS} -m
_move_runoff -c ${RUNOFF} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_command(TARGET pyroms_toolbox
POST_BUILD
COMMAND cp _average.so ${Project_SOURCE_DIR}/pyroms_toolbox/pyroms_toolbox/_average${CMAKE_SHARED_LIBRARY_SUFFIX}
COMMAND cp _move_runoff.so
${PROJECT_SOURCE_DIR}/pyroms_toolbox/pyroms_toolbox/_move_runoff${CMAKE_SHARED_LIBRARY_SUFFIX}
COMMAND cp -r ${Project_SOURCE_DIR}/pyroms_toolbox/pyroms_toolbox/ ${LIBRARY_OUTPUT_PATH}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
