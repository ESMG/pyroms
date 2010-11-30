CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(SRC
broyden.c
delaunay.c
geom.c
gridgen.c
hash.c
issimplepoly.c
istack.c
ode.c
swcr.c
vertlist.c
zode.c
)

set(LIBSRC
broyden.c
delaunay.c
geom.c
gridgen.c
hash.c
issimplepoly.c
istack.c
ode.c
swcr.c
vertlist.c
zode.c
c99-min/clog.c
c99-min/cmplx.c
c99-min/stubs.c
)

set(CFLAGS "-g -O2 -Wall -m64 -fPIC")
set(CFLAGS_TRIANGLE "-c -DTRILIBRARY -g -O2 -Wall -w -ffloat-store")
set(CFLAGS_GRIDGEN "-DGRIDGEN_STANDALONE ${CFLAGS}")

include_directories(.)
set_source_files_properties(triangle.c PROPERTIES COMPILE_FLAGS ${CFLAGS_TRIANGLE})
set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS ${CFLAGS})
add_executable(gridgen_executable main.c ${SRC} triangle.c)
set_target_properties(gridgen_executable PROPERTIES COMPILE_FLAGS ${CFLAGS_GRIDGEN})
target_link_libraries(gridgen_executable 
m 
"${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}gu${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

add_library(gridgen
SHARED
${LIBSRC} triangle.c
)
target_link_libraries(gridgen
m
"${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}gu${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

