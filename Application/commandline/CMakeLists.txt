
# add xThermal source and library
set(xThermal_DIR "${CMAKE_INSTALL_PREFIX}")
#include_directories("${xThermal_DIR}/include" "${xThermal_DIR}/include/PROST" "${xThermal_DIR}/include/CoolProp")
include_directories(${CMAKE_SOURCE_DIR}/Library ${CMAKE_SOURCE_DIR}/Library/thermo ${CMAKE_SOURCE_DIR}/Library/AMR_LUT ${CMAKE_SOURCE_DIR}/Library/H2ONaCl ${CMAKE_SOURCE_DIR}/Library/H2O ${CMAKE_SOURCE_DIR}/Library/NaCl)
#link_directories("${xThermal_DIR}/STATIC")
link_directories(${CMAKE_BINARY_DIR}/STATIC)
#link_directories(${CMAKE_BINARY_DIR}/SHARED)
aux_source_directory(. SRC_MAIN)
set( Srcs_bash ./SWEOSbash.cpp ./getopt_arguments.cpp)
set( Hdrs_bash ./SWEOSbash.h ./getopt_arguments.h)

# add Executable target
if(UNIX AND NOT APPLE)
    # fix static link issue for linux system
    set(LINK_FLAGS "-Wl,--no-as-needed -ldl")
endif()
add_executable(${PROJECT_NAME}_cmd  ${Srcs_bash}  main.cpp)
if (EXISTS ${LIBomp_static})
    target_link_libraries(${PROJECT_NAME}_cmd ${LIB_NAME} ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp_static})
else()
    target_link_libraries(${PROJECT_NAME}_cmd ${LIB_NAME} ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp})
endif ()

set_target_properties(${PROJECT_NAME}_cmd PROPERTIES OUTPUT_NAME ${PROJECT_NAME} RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} LINK_FLAGS "${LINK_FLAGS}")
install (TARGETS ${PROJECT_NAME}_cmd DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
