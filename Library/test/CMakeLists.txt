
list(APPEND TEST_ITEMS test_gsl test_omp test_PROST)
foreach(TEST_ITEM IN LISTS TEST_ITEMS)
  add_executable(${TEST_ITEM} ${TEST_ITEM}.cpp )
  target_link_libraries(${TEST_ITEM} ${GSL_LIBRARIES_STATIC} ${LIBomp} ${LIB_NAME}_STATIC)
  add_test(NAME ${TEST_ITEM} COMMAND ${TEST_ITEM})
  set_target_properties(${TEST_ITEM} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endforeach()