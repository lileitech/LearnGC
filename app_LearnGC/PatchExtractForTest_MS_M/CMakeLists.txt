#
# Add the executable
#

ADD_EXECUTABLE(PatchExtractForTest_MS_M PatchExtract_test.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h STLreader.cpp STLreader.h)
TARGET_LINK_LIBRARIES(PatchExtractForTest_MS_M zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh PatchExtractForTest_MS_M)
