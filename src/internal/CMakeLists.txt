#
# Libraries
#

#
# Binaries
#

# Triangle
ADD_EXECUTABLE(triangle triangle.c)
TARGET_LINK_LIBRARIES(triangle)
IF(UNIX)
  TARGET_LINK_LIBRARIES(triangle -lm)
ENDIF(UNIX)

SUBDIRS(io viewers)
