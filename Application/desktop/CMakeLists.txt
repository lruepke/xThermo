set(PROGRAM_NAME ${PROJECT_NAME}_desktop)
# add xThermal source and library
set(xThermal_DIR "${CMAKE_INSTALL_PREFIX}")
#include_directories("${xThermal_DIR}/include" "${xThermal_DIR}/include/PROST" "${xThermal_DIR}/include/CoolProp")
include_directories(${CMAKE_SOURCE_DIR}/Library ${CMAKE_SOURCE_DIR}/Library/thermo ${CMAKE_SOURCE_DIR}/Library/AMR_LUT ${CMAKE_SOURCE_DIR}/Library/H2ONaCl ${CMAKE_SOURCE_DIR}/Library/H2O ${CMAKE_SOURCE_DIR}/Library/NaCl)
#link_directories("${xThermal_DIR}/STATIC")
link_directories(${CMAKE_BINARY_DIR}/STATIC)
#link_directories(${CMAKE_BINARY_DIR}/SHARED)

macro(opt OPTION HELP VALUE)
    option(ENABLE_${OPTION} ${HELP} ${VALUE})
    set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro()
macro(set_config_option VARNAME STRING)
    set(${VARNAME} TRUE)
    list(APPEND CONFIG_OPTIONS ${STRING})
    message(STATUS "Found " ${STRING})
endmacro()

# global configuration
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)
# packager name
if(NOT xThermal_PACKAGER)
    execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT PACKAGER)
        set(PACKAGER "unknown")
    endif()
    string(REPLACE "\\" " " PACKAGER ${PACKAGER})
    set(xThermal_PACKAGER "${PACKAGER}")
    message(STATUS "Packager: " ${xThermal_PACKAGER})
endif()
# host system information
if(APPLE)
    set(xThermal_OS "MacOSX")
elseif(CYGWIN OR MSYS)
    # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
    # as a pure Windows build
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
            CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(CXX_COMPILER_MACHINE MATCHES "mingw")
            set(xThermal_OS "Windows")
            set(WIN32 1)
            add_definitions(-DWIN32)
        endif()
    endif()
else()
    set(xThermal_OS "${CMAKE_SYSTEM_NAME}")
endif()
message(STATUS "OS: " ${xThermal_OS})

# set icon
if(APPLE)
    set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Icons/logo.icns)
else()
    set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icon.rc)
endif()

# VTK and QT
set(Qt5_DIR=$ENV{Qt5_DIR})
message(STATUS "Qt5_DIR: " $ENV{Qt5_DIR})
if($ENV{Qt5_DIR})
    set(Qt5_DIR=$ENV{Qt5_DIR})
    message(STATUS "Qt5_DIR FOUND: " ${Qt5_DIR})
endif()
set(VTK_DIR=$ENV{VTK_DIR})
message(STATUS "VTK_DIR: " $ENV{VTK_DIR})
if($ENV{VTK_DIR})
    set(VTK_DIR=$ENV{VTK_DIR})
    message(STATUS "VTK_DIR FOUND: " ${VTK_DIR})
endif()
find_package(VTK COMPONENTS
        CommonCore
        FiltersCore
        InfovisCore
        InteractionStyle
        ViewsQt
        ChartsCore
        CommonDataModel
        RenderingContext2D
        RenderingContextOpenGL2
        RenderingCore
        RenderingFreeType
        RenderingGL2PSOpenGL2
        RenderingOpenGL2
        ViewsContext2D
        IOLegacy
        RenderingAnnotation
        )
if(VTK_FOUND)
    message(STATUS "VTK version: " "${VTK_VERSION}")
    message(STATUS "VTK QT version: " "${VTK_QT_VERSION}")
    # include(${VTK_USE_FILE}) #is no longer used starting with 8.90
    if("${VTK_VERSION}" STREQUAL "")
        message(FATAL_ERROR "Could not find VTK")
    else()
        set(xThermal_VTK_VERSION ${VTK_VERSION})
    endif()
    if("${VTK_QT_VERSION}" STREQUAL "")
        message(FATAL_ERROR "VTK was not built with Qt")
    else()
        set(xThermal_QT_VERSION ${VTK_QT_VERSION})
    endif()
else()
    message(FATAL_ERROR "VTK is not found")
endif()
#====================================================

# Use the include path and library for Qt that is used by VTK.
include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# GUI source files and resources
set( Srcs main.cxx  MainWindow.cxx qroundprogressbar.cpp  Meter.cpp)
set( Hdrs MainWindow.h qroundprogressbar.h  Meter.h)
set( MOC_Hdrs ui_MainWindow.h )
set( UIs MainWindow.ui )
set( QRCs Icons/icons.qrc)
set( Srcs_bash ../commandline/SWEOSbash.cpp ../commandline/getopt_arguments.cpp ${CMAKE_SOURCE_DIR}/Library/thermo/MultiProgressBar.cpp ${CMAKE_SOURCE_DIR}/Library/thermo/stdfunc.cpp)
set( Hdrs_bash ../commandline/SWEOSbash.h ${CMAKE_SOURCE_DIR}/Library/thermo/MultiProgressBar.h )
include_directories(../commandline) # commandline application path

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# The rest should just work (sure...)
# We have ui files, this will bring in the macro: qt5_wrap_ui
find_package(Qt5 COMPONENTS Widgets Core Concurrent LinguistTools REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QRC_Srcs ${QRCs} )

#languages
set(translation_zh "zh_CN")
set(TS_FILES languages/${translation_zh}.ts)
qt5_add_translation(QON_QM_FILES ${TS_FILES})

source_group("Resources" FILES
        ${UIs}
        ${QRCs}
        ${EXE_ICON} # Not present
        )
source_group("Generated" FILES
        ${UI_Srcs}
        ${MOC_Srcs}
        ${QRC_Srcs}
        )

# Add program
# if want to generage mac os app, add MACOSX_BUNDLE follow ${PROGRAM_NAME}
add_executable(
        ${PROGRAM_NAME}
        WIN32
        ${Srcs}
        ${Hdrs}
        ${Srcs_bash}
        ${Hdrs_bash}
        ${UI_Srcs}
        ${MOC_Hdrs}
        ${QRC_Srcs}
        ${myApp_ICON}
        ${QON_QM_FILES}
)
# test: disable macos bundle, cancle it after debug
# set(ENABLE_OS_SPECIFIC_INSTALL ON)
if(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
    message("Configure for macOS install")
    set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE ON)
    set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
    set_source_files_properties(${QON_QM_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources/translations")
    set(MACOSX_BUNDLE_ICON_FILE logo.icns)
    set_source_files_properties(${myApp_ICON} PROPERTIES
            MACOSX_PACKAGE_LOCATION "Resources")
endif()
if(UNIX AND NOT APPLE)
    # fix static link issue for linux system
    set(LINK_FLAGS "-Wl,--no-as-needed -ldl")
endif()
set_target_properties(${PROGRAM_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} LINK_FLAGS "${LINK_FLAGS}")
target_link_libraries(${PROGRAM_NAME} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Concurrent ${SWEOS_LIB_NAME} ${LIBFreeSteam} ${LIBomp})
target_link_libraries(${PROGRAM_NAME} ${VTK_LIBRARIES} )
vtk_module_autoinit(
        TARGETS ${PROGRAM_NAME}
        MODULES VTK::ChartsCore
        VTK::InfovisCore
        VTK::ViewsContext2D
        VTK::RenderingContextOpenGL2)
# link xThermal related libraries
#if (EXISTS ${LIBomp_static})
if(USE_OMP)
    target_link_libraries(${PROGRAM_NAME} ${LIB_NAME} ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp})
else()
    target_link_libraries(${PROGRAM_NAME} ${LIB_NAME} ${GSL_LIBRARIES_STATIC} ${CoolProp_LIBRARIES_STATIC} ${LIBomp_static})
endif ()


