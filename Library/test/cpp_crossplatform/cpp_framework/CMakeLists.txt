
# ==================================================================
#                generate cpp library framework
# ------------------------------------------------------------------
add_library(${PROJECT_NAME} SHARED example.cpp example.h)
set_target_properties(${PROJECT_NAME} PROPERTIES 
                    FRAMEWORK TRUE
                    FRAMEWORK_VERSION A
                    MACOSX_FRAMEWORK_IDENTIFIER ${FRAMEWORK_BUNDLE_IDENTIFIER}
                    MACOSX_FRAMEWORK_INFO_PLIST ${CMAKE_CURRENT_LIST_DIR}/framework.plist.in
                    # "current version" in semantic format in Mach-O binary file
                    VERSION 1.0.0
                    # "compatibility version" in semantic format in Mach-O binary file
                    SOVERSION 1.0.0
                    PUBLIC_HEADER "${CMAKE_CURRENT_LIST_DIR}/example.h"
                    XCODE_ATTRIBUTE_IPHONEOS_DEPLOYMENT_TARGET ${DEPLOYMENT_TARGET}
                    XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY ${CODE_SIGN_IDENTITY}
                    XCODE_ATTRIBUTE_DEVELOPMENT_TEAM ${DEVELOPMENT_TEAM_ID}
                    XCODE_ATTRIBUTE_TARGETED_DEVICE_FAMILY ${DEVICE_FAMILY}
                    XCODE_ATTRIBUTE_SKIP_INSTALL "YES"
                      )
install (TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
# ==================================================================
