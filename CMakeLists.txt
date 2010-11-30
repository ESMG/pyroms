CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)

add_subdirectory(external)
add_subdirectory(pyroms_toolbox)

add_custom_target(pyroms
COMMAND python ${Project_SOURCE_DIR}/pyinstall.py ${LIBRARY_OUTPUT_PATH}/../ 
WORKING_DIRECTORY ${Project_SOURCE_DIR}
)

