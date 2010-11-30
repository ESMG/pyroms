CMAKE_MINIMUM_REQUIRED(VERSION 2.4.6)

project(scrip Fortran)

if (CMAKE_Fortran_COMPILER MATCHES "pgf95$")
	set (F2PY_COMPILER "pg")
	set (CMAKE_FORTRAN_FLAGS "-g -r8 -Mnofree -O2 -fPIC")
elseif (CMAKE_Fortran_COMPILER MATCHES "pgf90$")
	set (F2PY_COMPILER "pg")
	set (CMAKE_FORTRAN_FLAGS "-g -r8 -Mnofree -O2 -fPIC")
elseif (CMAKE_Fortran_COMPILER MATCHES "gfortran$")
        set (F2PY_COMPILER "gnu95")
	set (CMAKE_FORTRAN_FLAGS "-g -fdefault-real-8 -ffixed-form -O2 -fPIC")
endif(CMAKE_Fortran_COMPILER MATCHES "pgf95$")

set(OBJF2PY
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/kinds_mod.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/constants.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/iounits.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/netcdf.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/grids.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_vars.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_distwgt.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_conserv.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_bilinear.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_bicubic.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip_test.dir/remap_read.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip_test.dir/remap.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/timers.o
${Project_BINARY_DIR}/external/scrip/source/CMakeFiles/scrip.dir/remap_write.o
)

set(OBJSET
kinds_mod.f
constants.f
iounits.f
netcdf.f
grids.f
remap_vars.f
remap_distwgt.f
remap_conserv.f
remap_bilinear.f
remap_bicubic.f
timers.f
remap_write.f
scrip.f
)

set(OBJTEST
kinds_mod.f
constants.f
iounits.f
netcdf.f
grids.f
timers.f
remap_vars.f
remap_read.f
remap.f
scrip_test.f
)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(Netcdf REQUIRED)
find_package(Netinc REQUIRED)

message("Netcdf Shared Library Found: ${Netcdf_FOUND}")
message("Netcdff Shared Library Found: ${Netcdff_FOUND}")
message("HDF5 Shared Library Found: ${Hdf5_FOUND}")
message("HDF5_HL Shared Library Found: ${Hdf5hl_FOUND}")

set_source_files_properties(${OBJSET} COMPILE_FLAGS ${CMAKE_FORTRAN_FLAGS})
set_source_files_properties(${OBJTEST} COMPILE_FLAGS ${CMAKE_FORTRAN_FLAGS})

include_directories(${Netcdf_INC})
add_executable(scrip
${OBJSET}
)

if ("${Netcdff_FOUND}" MATCHES "TRUE")
   if ("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
      target_link_libraries(scrip ${Netcdff_LIBRARY} ${Hdf5_LIBRARY} ${Hdf5hl_LIBRARY})
   else()
      target_link_libraries(scrip ${Netcdff_LIBRARY})
   endif("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
else()
   # This assumes that you are using a version of the NetCDF
   # shared library before version 4.0. This was before the
   # split of the shared library into the netcdf & netcdff
   # libraries. The Fortran functionality is needed for
   # scrip to build, and thus you must have the functions
   # available in either netcdf.so.3 or in netcdff.so
   target_link_libraries(scrip ${Netcdf_LIBRARY})
endif("${Netcdff_FOUND}" MATCHES "TRUE")

add_executable(scrip_test
${OBJTEST}
)
if ("${Netcdff_FOUND}" MATCHES "TRUE")
   if ("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
      target_link_libraries(scrip_test ${Netcdff_LIBRARY} ${Hdf5_LIBRARY} ${Hdf5hl_LIBRARY})
   else()
      target_link_libraries(scrip_test ${Netcdff_LIBRARY})
   endif("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
else()
   # Again, this assumes that you have an old version of the
   # Netcdf shared library object. It needs to be before version
   # 4.0 since this is when the Fortran functionality was removed
   # from the netcdf shared library and placed into the netcdff
   # shared library.
   target_link_libraries(scrip_test ${Netcdf_LIBRARY})
endif("${Netcdff_FOUND}" MATCHES "TRUE")

#add_custom_command(TARGET scrip_test 
#POST_BUILD	
#COMMAND cp ${CMAKE_SOURCE_DIR}/external/scrip/source/CMakeFiles/scrip.dir/*.o ${CMAKE_SOURCE_DIR}/external/scrip/source/
#COMMAND cp
#${CMAKE_SOURCE_DIR}/external/scrip/source/CMakeFiles/scrip_test.dir/*.o ${CMAKE_SOURCE_DIR}/external/scrip/source/
#)

set(pyscrip "${Project_SOURCE_DIR}/external/scrip/source/pyscrip.f90")

string(REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}netcdf${CMAKE_SHARED_LIBRARY_SUFFIX}" "" Netcdf_DIR ${Netcdf_LIBRARY})

if ("${Netcdff_FOUND}" MATCHES "TRUE")
   string(REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}netcdff${CMAKE_SHARED_LIBRARY_SUFFIX}" "" Netcdff_DIR ${Netcdff_LIBRARY})
   if ("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
       string(REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}hdf5${CMAKE_SHARED_LIBRARY_SUFFIX}" "" Hdf5_DIR ${Hdf5_LIBRARY})

       string(REPLACE "${CMAKE_SHARED_LIBRARY_PREFIX}hdf5_hl${CMAKE_SHARED_LIBRARY_SUFFIX}" "" Hdf5hl_DIR ${Hdf5hl_LIBRARY})

       add_custom_target(scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
       COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${CMAKE_FORTRAN_FLAGS} -L${Netcdf_DIR} -lnetcdf -L${Netcdff_DIR} -lnetcdff -L${Hdf5_DIR} -lhdf5 -L${Hdf5hl_DIR} -lhdf5_hl -m scrip -c ${OBJF2PY} ${pyscrip} 
       )
   else()
       add_custom_target(scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
       COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${CMAKE_FORTRAN_FLAGS} -L${Netcdf_DIR} -lnetcdf -L${Netcdff_DIR} -lnetcdff -m scrip -c ${OBJF2PY} ${pyscrip} 
       )
   endif("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
else()
   if ("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
      add_custom_target(scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
      COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${CMAKE_FORTRAN_FLAGS} -L${Netcdf_DIR} -lnetcdf -L${Hdf5_DIR} -hdf5 -L${Hdf5hl_DIR} -lhdf5_hl -m scrip -c ${OBJF2PY} ${pyscrip})
   else()
      add_custom_target(scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
      COMMAND f2py --fcompiler=${F2PY_COMPILER} --f90flags=${CMAKE_FORTRAN_FLAGS} -L${Netcdf_DIR} -lnetcdf -m scrip -c ${OBJF2PY} ${pyscrip}
      )
   endif("${Hdf5_FOUND}" MATCHES "TRUE" AND "${Hdf5hl_FOUND}" MATCHES "TRUE")
endif("${Netcdff_FOUND}" MATCHES "TRUE") 

add_dependencies(scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
scrip_test
)

if("scrip.so" MATCHES "scrip${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_custom_command(TARGET scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
POST_BUILD
COMMAND cp ${Project_SOURCE_DIR}/external/scrip/source/scrip${CMAKE_SHARED_LIBRARY_SUFFIX} ${LIBRARY_OUTPUT_PATH}
)
ELSE("scrip.so" MATCHES "scrip${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_custom_command(TARGET scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
POST_BUILD
COMMAND mv ${Project_SOURCE_DIR}/external/scrip/source/scrip.so ${Project_SOURCE_DIR}/external/scrip/source/scrip${CMAKE_SHARED_LIBRARY_SUFFIX}
COMMAND cp ${Project_SOURCE_DIR}/external/scrip/source/scrip${CMAKE_SHARED_LIBRARY_SUFFIX} ${LIBRARY_OUTPUT_PATH}
)
ENDIF("scrip.so" MATCHES "scrip${CMAKE_SHARED_LIBRARY_SUFFIX}")


