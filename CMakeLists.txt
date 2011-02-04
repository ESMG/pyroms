CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)

if( NOT APPLE )
# check for 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
   # If the size of a void pointer is 4 bytes, this is a 32 bit machine
   # Set have 64 bit to zero
   set(HAVE_64_BIT 0)
else( CMAKE_SIZEOF_VOID_P EQUAL 4)
   # If the size of a void pointer is > 4 bytes (8 bytes), this is a 64 bit
   # machine. Set have 64 bit to one.
   set(HAVE_64_BIT 1)
endif( CMAKE_SIZEOF_VOID_P EQUAL 4)
endif( NOT APPLE )

#configure_file( ${Project_SOURCE_DIR}/config.h.in
#${EXECUTABLE_OUTPUT_PATH}/config.h
#)
#include_directories(${EXECUTABLE_OUTPUT_PATH})

add_subdirectory(external)
add_subdirectory(pyroms_toolbox)

set(OUTPUT ${LIBRARY_OUTPUT_PATH}/../)
string(REGEX MATCH "/lib/../$" IS_LIB ${OUTPUT})
if(${IS_LIB} MATCHES "/lib/../")
add_custom_target(pyroms
COMMAND python ${Project_SOURCE_DIR}/pyinstall.py ${LIBRARY_OUTPUT_PATH}/../ 
WORKING_DIRECTORY ${Project_SOURCE_DIR}
)
else(${IS_LIB} MATCHES "/lib/../")
add_custom_target(pyroms
COMMAND python ${Project_SOURCE_DIR}/pyinstall.py ${LIBRARY_OUTPUT_PATH}/
WORKING_DIRECTORY ${Project_SOURCE_DIR}
)
endif(${IS_LIB} MATCHES "/lib/../")

