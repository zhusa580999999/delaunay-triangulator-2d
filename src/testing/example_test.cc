// An example implementation of a test suite.

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class ExampleTest : public CPPUNIT_NS::TestCase {
  CPPUNIT_TEST_SUITE(ExampleTest);
  CPPUNIT_TEST(testAddition);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp() {}
  void tearDown() {} 

 protected:
  void testAddition() {
   CPPUNIT_ASSERT(5 == 2 + 3);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ExampleTest);
