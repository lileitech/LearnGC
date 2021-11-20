#
# Add the executable
#

ADD_EXECUTABLE(PatchExtractForTest_MS PatchExtract_test.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(PatchExtractForTest_MS zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh PatchExtractForTest_MS)
