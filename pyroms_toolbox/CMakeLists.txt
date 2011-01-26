CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)

add_subdirectory(pyroms_toolbox)

add_custom_command(TARGET pyroms_toolbox
POST_BUILD
COMMAND cp -r ${Project_SOURCE_DIR}/pyroms_toolbox/pyroms_toolbox ${LIBRARY_OUTPUT_PATH}
)
