CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(LIBSRC
gridmap.c
gridnodes.c
gucommon.c
poly.c
)

set(LIBHEADERS
gridmap.h
gridnodes.h
gucommon.h
poly.h
)

set(SRC
getbound.c
getnodes.c
insertgrid.c
setbathy.c
subgrid.c
xy2ij.c
gridaverager.c
gridbathy.c
)

set(OBJ
getbound
getnodes
insertgrid
setbathy
subgrid
xy2ij
)


if(HAVE_64_BIT)
   set(CFLAGS "-g -O2 -Wall -pedantic -m64 -fPIC")
else(HAVE_64_BIT)
   set(CFLAGS "-g -O2 -Wall -pedantic -fPIC")
endif(HAVE_64_BIT)

set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "${CFLAGS}")

add_executable(getbound
getbound.c 
)

add_executable(getnodes
getnodes.c
)

add_executable(insertgrid
insertgrid.c
)

add_executable(setbathy
setbathy.c
)

add_executable(subgrid
subgrid.c
)

add_executable(xy2ij
xy2ij.c
)

set(BATHYSRC
gridaverager.c
gridbathy.c
)

include_directories("${Project_SOURCE_DIR}/external/nn")
include_directories("${Project_SOURCE_DIR}/external/csa")

message("${Project_SOURCE_DIR}")

add_executable(gridbathy
${BATHYSRC}
)

add_library(gu
SHARED
${LIBSRC}
${LIBHEADERS}
)

target_link_libraries(getbound gu m)
target_link_libraries(getnodes gu m)
target_link_libraries(insertgrid gu m)
target_link_libraries(setbathy gu m)
target_link_libraries(subgrid gu m)
target_link_libraries(xy2ij gu m)
target_link_libraries(gridbathy 
gu 
m 
"${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}nn${CMAKE_SHARED_LIBRARY_SUFFIX}"
"${LIBRARY_OUTPUT_PATH}/${CMAKE_STATIC_LIBRARY_PREFIX}csa${CMAKE_STATIC_LIBRARY_SUFFIX}"
) 

