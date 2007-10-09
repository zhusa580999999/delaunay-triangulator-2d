IF(X11_FOUND)
  ADD_EXECUTABLE(showme
                 showme.c)
  TARGET_LINK_LIBRARIES(showme
                        ${X11_LIBRARIES})
ENDIF(X11_FOUND)
