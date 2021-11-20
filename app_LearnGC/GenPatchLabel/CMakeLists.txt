#
# Add the executable
#

ADD_EXECUTABLE(GenPatchLabel Generate_PatchLabel.cpp STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(GenPatchLabel zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh GenPatchLabel)
