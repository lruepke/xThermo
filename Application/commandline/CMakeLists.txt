# Command-line tool "xThermal"
add_executable(xThermal_cli main.cpp SWEOSbash.cpp getopt_arguments.cpp)
target_link_libraries(xThermal_cli PRIVATE xThermal)
target_compile_options(xThermal_cli PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-writable-strings>)
# the executable sits in its own directory: in the build tree, "xThermal" next to the Python package directory
# "python/xThermal" would collide on case-insensitive file systems
set_target_properties(xThermal_cli PROPERTIES OUTPUT_NAME xThermal RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
install(TARGETS xThermal_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
