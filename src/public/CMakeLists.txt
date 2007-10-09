#
# Libraries
#

# Triangle
ADD_LIBRARY(triangle ../internal/triangle.c)
SET_TARGET_PROPERTIES(triangle PROPERTIES COMPILE_FLAGS "-DTRILIBRARY")
IF(UNIX)
  TARGET_LINK_LIBRARIES(triangle -lm)
ENDIF(UNIX)

#
# Binaries
#

# Utility to show how the Triangle library is to be used
ADD_EXECUTABLE(tricall tricall.c)
TARGET_LINK_LIBRARIES(tricall triangle)
IF(UNIX)
  TARGET_LINK_LIBRARIES(triangle -lm)
ENDIF(UNIX)

