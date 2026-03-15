if(NOT Build_API_Julia)
    return()
endif(NOT Build_API_Julia)

# ---------------------------- Find Julia ---------------------------------
find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
# It seems that the Julia only support dynamic lib wrapper.
# -------------------------------------------------------------------------
set(INSTALL_PATH_API_Julia ${CMAKE_INSTALL_PREFIX}/API/julia/${PROJECT_NAME})

function(Build_Julia_API Module SC_FILES libs_depend)
    message(STATUS "Build Julia module ${Module}")
    add_library(${Module}_jl SHARED ${SC_FILES} ${Module}.cxx) #cxx is exclusively used as Julia warpper extension
    # set_target_properties("${LIB_NAME}_${buildType}" PROPERTIES OUTPUT_NAME ${LIB_NAME})
    target_link_libraries(${Module}_jl JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl ${libs_depend})
    # install
    install(TARGETS ${Module}_jl DESTINATION ${INSTALL_PATH_API_Julia}/${Module})
    INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/${Module}/ DESTINATION ${INSTALL_PATH_API_Julia}
      FILES_MATCHING PATTERN "*.jl"
      PATTERN "IAPWS_Others" EXCLUDE
      )
endfunction(Build_Julia_API)

