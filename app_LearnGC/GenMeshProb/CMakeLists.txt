#
# Add the executable
#

ADD_EXECUTABLE(GenMeshProb Generate_MeshProb.cpp STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(GenMeshProb zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh GenMeshProb)
