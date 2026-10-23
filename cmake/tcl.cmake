if(NOT Build_API_tcl)
    return()
endif(NOT Build_API_tcl)

# ---------------------------- Find Julia ---------------------------------
FIND_PACKAGE(Tcl)  
if(TCL_FOUND)
    message("TCL include found: " ${TCL_INCLUDE_PATH} )
    message("TCL libs found: " ${TCL_LIBRARY} )
    INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH})
    SET(CMAKE_SWIG_FLAGS "")
endif(TCL_FOUND)
# -------------------------------------------------------------------------
set(INSTALL_PATH_API_tcl ${CMAKE_INSTALL_PREFIX}/API/tcl/${PROJECT_NAME})

function(Build_tcl_API Module API_SRC)
    swig_add_library(tcl${Module}
    LANGUAGE tcl
    SOURCES ${Module}.i ${API_SRC})
    SWIG_LINK_LIBRARIES(tcl${Module} ${TCL_LIBRARY})
    # install
    set_target_properties( tcl${Module} PROPERTIES OUTPUT_NAME "${Module}" SUFFIX "")
    install(TARGETS tcl${Module} DESTINATION ${INSTALL_PATH_API_tcl} "${Module}")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/test_${Module}.tcl DESTINATION ${INSTALL_PATH_API_tcl})
endfunction(Build_tcl_API )
