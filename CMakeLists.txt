project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-g -Wall -Wno-int-in-bool-context")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(src_oo "src/oo_version")

set(sources src/main.cpp src/tools_fp.cpp src/FusionEKF_fp.cpp  src/tools_fp.h src/FusionEKF_fp.h)

set(sources_test src/test_fp.cpp src/tools_fp.cpp src/FusionEKF_fp.cpp  src/tools_fp.h src/FusionEKF_fp.h)

set(sources_oo ${src_oo}/main.cpp ${src_oo}/tools.cpp ${src_oo}/FusionEKF.cpp ${src_oo}/kalman_filter.cpp ${src_oo}/tools.h src/FusionEKF.h ${src_oo}/kalman_filter.h)
set(sources_test_oo ${src_oo}/test.cpp ${src_oo}/tools.cpp ${src_oo}/FusionEKF.cpp ${src_oo}/kalman_filter.cpp ${src_oo}/tools.h ${src_oo}/FusionEKF.h ${src_oo}/kalman_filter.h)

message( "System name = ${CMAKE_SYSTEM_NAME}" )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

add_executable(ExtendedKF ${sources})
add_executable(my_test ${sources_test})

#add_executable(my_test_oo ${sources_test_oo})
target_link_libraries(ExtendedKF z ssl uv uWS)
target_link_libraries(my_test z ssl uv uWS)
