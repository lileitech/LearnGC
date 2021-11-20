#
# Add the executable
#

ADD_EXECUTABLE(ScarSegByThreshold main.cpp STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(ScarSegByThreshold  zxhImages zxhCommon zxhImageModelings ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh ScarSegByThreshold)
