project (libColorConvert CXX)

SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -g -ggdb -O3 -fopenmp")

FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_Include_Dirs})

add_library(ColorConvert colorconvert.cpp)

add_executable(testColorConvert testcolorconvert.cpp)
target_link_libraries(testColorConvert ColorConvert)

target_link_libraries(ColorConvert ${OpenCV_LIBS})

find_package(GTest REQUIRED)
target_link_libraries(testColorConvert ${GTEST_BOTH_LIBRARIES})
