#
# Add the executable
#

ADD_EXECUTABLE(GenLAlabel_zxh GenerateLALabel_zxh.cpp)
TARGET_LINK_LIBRARIES(GenLAlabel_zxh zxhImageModelings zxhImages zxhCommon ${NIFTI_NIFTILIB_NAME} )
 
INSTALL_TARGETS(/bin/zxh GenLAlabel_zxh)
