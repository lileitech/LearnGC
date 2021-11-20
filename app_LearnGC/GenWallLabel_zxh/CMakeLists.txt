#
# Add the executable
#

ADD_EXECUTABLE(GenWallLabel_zxh GenerateWallLabel_zxh.cpp STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(GenWallLabel_zxh zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh GenWallLabel_zxh)
