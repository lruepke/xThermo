if(USE_COOLPROP)
    add_compile_definitions("USE_COOLPROP=1" )
    set(CoolProp_DIR "${CMAKE_SOURCE_DIR}/ThirdParties/install/CoolProp_${TARGET_SUFFIX}" CACHE FILEPATH "CoolProp path contains include, lib")
    if("${CoolProp_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}/ThirdParties/install/CoolProp_${TARGET_SUFFIX}")
        message(STATUS "CoolProp_DIR is not changed")
    else()
        unset(CoolProp_DIR CACHE)
        set(CoolProp_DIR "${CMAKE_SOURCE_DIR}/ThirdParties/install/CoolProp_${TARGET_SUFFIX}" CACHE FILEPATH "CoolProp path contains include, lib")
    endif()
    if(EXISTS ${CoolProp_DIR})
        set(CoolProp_INCLUDE_DIR "${CoolProp_DIR}/include")
        INCLUDE_DIRECTORIES(${CoolProp_INCLUDE_DIR})
        if(WIN32)
            if(BUILD_SHARED)
                add_compile_definitions("COOLPROP_DLL=1")
                set(CoolProp_LIBRARIES_SHARED ${CoolProp_DIR}/shared/CoolProp.lib )
                set(CoolProp_LIBRARIES_STATIC ${CoolProp_DIR}/shared/CoolProp.lib)
            else()
                set(CoolProp_LIBRARIES_SHARED ${CoolProp_DIR}/shared/CoolProp.lib )
                set(CoolProp_LIBRARIES_STATIC ${CoolProp_DIR}/static/CoolProp.lib)
            endif()
        else()
            set(CoolProp_LIBRARY_DIR "${CoolProp_DIR}/lib")
            set(CoolProp_LIBRARIES_SHARED ${CoolProp_LIBRARY_DIR}/libCoolProp${CMAKE_SHARED_LIBRARY_SUFFIX} )
            set(CoolProp_LIBRARIES_STATIC ${CoolProp_LIBRARY_DIR}/libCoolProp.a)
            link_directories(${CoolProp_LIBRARY_DIR})
            FILE(GLOB libfiles "${CoolProp_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
        endif()
        INSTALL(DIRECTORY ${CoolProp_INCLUDE_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/CoolProp)
        INSTALL(FILES ${CoolProp_LIBRARIES_STATIC} DESTINATION ${CMAKE_INSTALL_PREFIX}/STATIC)
        INSTALL(FILES ${CoolProp_LIBRARIES_SHARED} DESTINATION ${CMAKE_INSTALL_PREFIX}/SHARED)
        message(STATUS "CoolProp include path found: " ${CoolProp_INCLUDE_DIR} )
        message(STATUS "CoolProp library path found: " ${CoolProp_LIBRARY_DIR})
        # copy coolprop to current binary path
        function(CopyCoolPropLib name_target)
            if(WIN32)
                message(STATUS "Cannot copy CoolProp.dll to the project binary path.")
            else()
                add_custom_command(
                    TARGET
                    ${name_target}
                    POST_BUILD COMMAND /bin/sh -c
                    \"COMMAND_DONE=0 \;
                    if ${CMAKE_COMMAND} -E copy_if_different
                        \\${CoolProp_LIBRARIES_SHARED}
                        \\${CMAKE_BINARY_DIR}
                        \&\>/dev/null \; then
                        COMMAND_DONE=1 \;
                    fi \;
                    if [ \\$$COMMAND_DONE -eq 0 ] \; then
                        echo Failed to copy the \\${CoolProp_LIBRARIES_SHARED} to binary directory \\${CMAKE_BINARY_DIR} \;
                        exit 1 \;
                    fi\"
                )
            endif()
        endfunction(CopyCoolPropLib)
    else()
        message(FATAL_ERROR "CoolProp_DIR = ${CoolProp_DIR} doesn't exist. Please use option -DCoolProp_DIR=xxx set a correct path which contains include and lib dir")
    endif()
    INCLUDE_DIRECTORIES(${CoolProp_INCLUDE_DIR})
endif()
