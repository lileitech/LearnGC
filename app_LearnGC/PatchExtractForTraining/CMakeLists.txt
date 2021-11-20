#
# Add the executable
#

ADD_EXECUTABLE(PatchExtractForTraining PatchExtract_train.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(PatchExtractForTraining zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh PatchExtractForTraining)
