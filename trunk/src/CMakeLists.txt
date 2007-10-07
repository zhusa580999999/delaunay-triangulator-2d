INCLUDE(${CMAKE_ROOT}/Modules/FindOpenGL.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/FindGLUT.cmake)

INCLUDE_DIRECTORIES(
  ${TRIANGLE_SOURCE_DIR}/inc
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${X11_INCLUDE_DIR})

ADD_LIBRARY(reader SHARED
  spreader_spa.cpp
  spreader_spb.cpp
  spreader_node.cpp
  spreader_raw.cpp
  spreader_raw_d.cpp
  spreader_tiles.cpp
  smreader_sma.cpp
  smreader_smb.cpp
  svreader_sva.cpp
  svreader_svb.cpp)

ADD_LIBRARY(writer SHARED
  spwriter_spa.cpp
  spwriter_spb.cpp
  spwriter_raw.cpp
  smwriter_sma.cpp
  smwriter_smb.cpp)

ADD_LIBRARY(cspglue SHARED 
  csp_glue.cpp)

TARGET_LINK_LIBRARIES(cspglue
  reader
  writer)

ADD_EXECUTABLE(triangle
  triangle.c
  triextra.c)
TARGET_LINK_LIBRARIES(triangle
  cspglue)
IF(UNIX)
  TARGET_LINK_LIBRARIES(triangle
    -lm)
ENDIF(UNIX)

IF(OPENGL_FOUND AND OPENGL_GLU_FOUND AND GLUT_FOUND)
  ADD_EXECUTABLE(spviewer
    delaunay3.cpp
    delaunay2.cpp
    predicates.cpp
    spdelaunay2d.cpp
    spdelaunay3d.cpp
    spreadscattered.cpp
    sscontainer3d.cpp
    sscontainer2d.cpp
    sp_viewer.cpp
    spconverter.cpp)
  TARGET_LINK_LIBRARIES(spviewer
    reader
    ${OPENGL_LIBRARIES}
    ${GLUT_glut_LIBRARY})
ENDIF(OPENGL_FOUND AND OPENGL_GLU_FOUND AND GLUT_FOUND)

IF(X11_FOUND)
  ADD_EXECUTABLE(showme
    showme.c)
  TARGET_LINK_LIBRARIES(showme
    ${X11_LIBRARIES})
ENDIF(X11_FOUND)

ADD_EXECUTABLE(sp2sp
  sp2sp.cpp)
TARGET_LINK_LIBRARIES(sp2sp
  reader
  writer)

FOREACH(DataSet ${DataSets})
  ADD_TEST(${DataSet}_clean rm ${DataDir}/${DataSet}.1.sma ${DataDir}/${DataSet}.1.spa)
  ADD_TEST(${DataSet}_finalize ${EXECUTABLE_OUTPUT_PATH}/triangle -f ${DataDir}/${DataSet}.spb)
  ADD_TEST(${DataSet}_check_header diff ${DataDir}/${DataSet}.1.sma ${TestDataDir}/${DataSet}.1.sma)
  ADD_TEST(${DataSet}_check_point_stream diff ${DataDir}/${DataSet}.1.spa ${TestDataDir}/${DataSet}.1.spa)
ENDFOREACH(DataSet)
