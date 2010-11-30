CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

set(LIBSRC
csa.c
svd.c
)

set(LIBHEADERS
csa.h
minell.h
nan.h
svd.h
version.h
)

set(OBJECTS
csa
minell
svd
)

set(SRC
csa.c
minell.c
svd.c
)

set(CFLAGS "-g -O2 -Wall -pedantic")
set(CFLAGS_BATHY "${CFLAGS} -DCSA_STANDALONE")

set_source_files_properties(${LIBSRC} PROPERTIES COMPILE_FLAGS ${CFLAGS})

add_executable(csabathy 
csa 
minell 
svd
)
set_target_properties(csabathy PROPERTIES COMPILE_FLAGS ${CFLAGS_BATHY})
target_link_libraries(csabathy m) 

add_library(csa
STATIC
csa
svd
)

target_link_libraries(csa m)
