# 一些关于苹果系统的设置

# set compiler flag
if (APPLE)
  add_definitions("-Wno-array-bounds -Wno-unused-command-line-argument") # turn off this array bound check warning at this moment, this warning accurs at where dim!=3
  set(CMAKE_OSX_DEPLOYMENT_TARGET 12.0) #used to get rid of the annoying warning: *** was built for newer macOS version (12.0) than being linked (11.3)
endif()