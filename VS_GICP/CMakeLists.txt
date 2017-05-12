cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(gicp_alg)

find_package(PCL 1.8 REQUIRED)

set(MY_INCLUDE_DIR "./")
set(MATLAB_INCLUDE_DIR "$ENV{MATLAB_ROOT}/extern/include")
set(MATLAB_LIB_DIR "$ENV{MATLAB_ROOT}/extern/lib/win64/microsoft")

include_directories(${PCL_INCLUDE_DIRS} ${MY_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR})
link_directories(${PCL_LIBRARY_DIRS} ${MATLAB_LIB_DIR})
add_definitions(${PCL_DEFINITIONS} /D "MATLAB_MEX_FILE")

add_library(gicp_alg SHARED gicp_main.cpp def_file.def)
set_target_properties(gicp_alg PROPERTIES SUFFIX ".mexw64")
target_link_libraries(gicp_alg ${PCL_LIBRARIES} "libmx.lib" "libmex.lib" "libmat.lib")
add_custom_command(TARGET gicp_alg POST_BUILD COMMAND xcopy $(TargetPath) $(ProjectDir)..\\.. /Y)