INCLUDE(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
# INCLUDE(${CMAKE_ROOT}/Modules/FindGLUT.cmake)   # for spviewer

INCLUDE_DIRECTORIES(
#  ${OPENGL_INCLUDE_DIR}   # for spviewer
#  ${GLUT_INCLUDE_DIR}     # for spviewer
  ${X11_INCLUDE_DIR})

IF(X11_FOUND)
  ADD_EXECUTABLE(showme
                 showme.c)
  TARGET_LINK_LIBRARIES(showme
                        ${X11_LIBRARIES})
ENDIF(X11_FOUND)
