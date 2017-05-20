project(GoICP)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(MY_INCLUDE_DIR "./")
set(MATLAB_INCLUDE_DIR "$ENV{MATLAB_ROOT}/extern/include")
set(MATLAB_LIB_DIR "$ENV{MATLAB_ROOT}/extern/lib/win64/microsoft")

include_directories(${MY_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR})

link_directories(${MATLAB_LIB_DIR})

add_definitions(/D "MATLAB_MEX_FILE")

add_library(GoICP SHARED
	mex_main.cpp
	jly_goicp.cpp
	jly_3ddt.cpp
	matrix.cpp
	ConfigMap.cpp
	StringTokenizer.cpp
	def_file.def
	)

set_target_properties(GoICP PROPERTIES SUFFIX ".mexw64")

target_link_libraries(GoICP "libmx.lib" "libmex.lib" "libmat.lib")

add_custom_command(TARGET GoICP POST_BUILD COMMAND xcopy $(TargetPath) $(ProjectDir)..\\.. /Y)
