INCLUDE(FindCppUnit.cmake)

#
# Libraries
#

# CppUnit + main().
#
# Link this target into the code containing your testing fixtures.
# See example_test for an example.
ADD_LIBRARY(testing_main testing_main.cc)
TARGET_LINK_LIBRARIES(testing_main ${CPPUNIT_LIBRARIES})

#
# Tests
#

# Example test suite
ADD_EXECUTABLE(example_test example_test.cc)
ADD_TEST(example_test ${EXECUTABLE_OUTPUT_PATH}/example_test)
TARGET_LINK_LIBRARIES(example_test testing_main)

