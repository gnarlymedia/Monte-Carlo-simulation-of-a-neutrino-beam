cmake_minimum_required(VERSION 3.6)
project(proj_3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES proj_3.c)

add_executable(proj_3 ${SOURCE_FILES})

find_library(M_LIB m X11)
find_library(CPGPLOT_LIB cpgplot pgplot png z)
target_link_libraries(proj_3 ${M_LIB})
target_link_libraries(proj_3 ${CPGPLOT_LIB})