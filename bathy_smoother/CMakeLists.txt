CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# This variable must be set to the location for binary files 
#set(BIN_LOCATION "/u1/uaf/torgerso/igloo/bin/")

# This variable must be set to the location of your library files location
#set(LIB_LOCATION "/u1/uaf/torgerso/igloo/lib/")

# This variable must be set to the location of your Python code in your PYTHONPATH
#set(PYTHON_PATH "/u1/uaf/torgerso/igloo/python/")

add_subdirectory(external)

add_custom_command(OUTPUT lpsolve
COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lp_solve/ccc
COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lp_solve/bin/ux64/lp_solve ${EXECUTABLE_OUTPUT_PATH}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lp_solve/
)

add_custom_command(OUTPUT lpsolve55
COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lpsolve55/ccc
COMMAND cp
"${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lpsolve55/bin/ux64/${CMAKE_STATIC_LIBRARY_PREFIX}lpsolve55${CMAKE_STATIC_LIBRARY_SUFFIX}" ${LIBRARY_OUTPUT_PATH}
COMMAND cp
"${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lpsolve55/bin/ux64/${CMAKE_SHARED_LIBRARY_PREFIX}lpsolve55${CMAKE_SHARED_LIBRARY_SUFFIX}" ${LIBRARY_OUTPUT_PATH}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/lpsolve55/
)

add_custom_target(pylibs
COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/extra/Python/pyhelp.py ${LIBRARY_OUTPUT_PATH}/
COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/extra/Python/setup.py install --prefix="${LIBRARY_OUTPUT_PATH}/"
DEPENDS lpsolve lpsolve55
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/lp_solve_5.5/extra/Python/
)

add_custom_target(bathy_smoother
COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/setup.py install
--prefix="${LIBRARY_OUTPUT_PATH}/"
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
add_dependencies(bathy_smoother pylibs)
