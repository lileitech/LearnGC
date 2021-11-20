#
# Add the executable
#

ADD_EXECUTABLE(GraphCut ComputeTime.cpp ComputeTime.h graph.cpp graph.h maxflow.cpp main.cpp segmenter.h STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(GraphCut  zxhImages zxhCommon zxhImageModelings ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh GraphCut)
