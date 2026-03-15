find_package (Git)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                  OUTPUT_VARIABLE GIT_TAG_LATEST
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
string(TIMESTAMP COMPILE_TIME %Y-%m-%d) #_%H%M%S
message(STATUS "The project latest version: " ${GIT_TAG_LATEST})
string(REGEX REPLACE "^[Vv]([0-9]+)\\..*" "\\1" VERSION_MAJOR "${GIT_TAG_LATEST}")
string(REGEX REPLACE "^[Vv][0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${GIT_TAG_LATEST}")
string(REGEX REPLACE "^[Vv][0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${GIT_TAG_LATEST}")

set(xThermal_MAJOR_VERSION ${VERSION_MAJOR})
set(xThermal_MINOR_VERSION ${VERSION_MINOR})
set(xThermal_PATCH_VERSION ${VERSION_PATCH})
set(xThermal_EXTRA_VERSION "")
set(xThermal_EXTRA_VERSION_TEXI "${xThermal_EXTRA_VERSION}")
if(NOT xThermal_RELEASE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(GIT_COMMIT_HASH)
    set(xThermal_EXTRA_VERSION "${xThermal_EXTRA_VERSION}-git-${GIT_COMMIT_HASH}")
  endif()
    set(xThermal_EXTRA_VERSION_TEXI "${xThermal_EXTRA_VERSION_TEXI} (development version)")
endif()
set(xThermal_VERSION "${xThermal_MAJOR_VERSION}.${xThermal_MINOR_VERSION}")
set(xThermal_VERSION "${xThermal_VERSION}.${xThermal_PATCH_VERSION}${xThermal_EXTRA_VERSION}")
set(xThermal_DATE ${COMPILE_TIME})
set(xThermal_SHORT_LICENSE "GNU General Public License")
message(STATUS "TARGET_SUFFIX: ${TARGET_SUFFIX}")