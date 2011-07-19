CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

set(LIBSRC 
delaunay.c 
hash.c 
istack.c 
lpi.c 
minell.c 
nnai.c 
nnpi.c 
nncommon.c 
nncommon-vulnerable.c 
preader.c
)

set(OBJECTS
delaunay
hash
istack
lpi
minell
nnai
nnpi
nncommon
)

set(LIBHEADERS
delaunay.h
hash.h
istack.h
minell.h
nan.h
nn.h
nn_internal.h
preader.h
version.h
)

set(CFLAGS "-g -O2 -Wall -pedantic")

set(CFLAGS_TRIANGLE "-c -DTRILIBRARY -g -O2 -Wall -pedantic -w -ffloat-store")

set(CFLAGS_NNBATHY "${CFLAGS}") 

set(CFLAGS_MINELL "-DME_STANDALONE ${CFLAGS}")

set(CFLAGS_VUL "${CFLAGS} -ffloat-store")

set_source_files_properties(${LIBSRC} PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_source_files_properties(triangle.c PROPERTIES COMPILE_FLAGS ${CFLAGS_TRIANGLE})
set_source_files_properties(nncommon-vulnerable.c PROPERTIES COMPILE_FLAGS ${CFLAGS_VUL})

include_directories(.)

add_library(nn
SHARED
triangle.c ${LIBSRC} 
)

message("${LIBRARY_OUTPUT_PATH}")
link_directories(.)

set_source_files_properties(minell.c PROPERTIES COMPILE_FLAGS ${CFLAGS_MINELL})
add_executable(minell minell.c)
target_link_libraries(minell 
m
)
add_dependencies(minell pyroms)

set_source_files_properties(nnbathy.c PROPERTIES COMPILE_FLAGS ${CFLAGS_NNBATHY})
add_executable(nnbathy nnbathy.c)
target_link_libraries(nnbathy
"${LIBRARY_OUTPUT_PATH}/${CMAKE_SHARED_LIBRARY_PREFIX}nn${CMAKE_SHARED_LIBRARY_SUFFIX}"
m
)
