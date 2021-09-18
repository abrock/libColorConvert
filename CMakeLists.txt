project (libColorConvert CXX)

SET(CMAKE_CXX_FLAGS " -Wall -pedantic -g -ggdb -O3 -fopenmp -Werror")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_Include_Dirs})

add_library(ColorConvert
    colorconvert.cpp
    colorflow.cpp
    )

add_executable(testColorConvert testcolorconvert.cpp)
target_link_libraries(testColorConvert ColorConvert)

add_executable(largest-circle largest-circle.cpp)
target_link_libraries(largest-circle ColorConvert)

target_link_libraries(ColorConvert ${OpenCV_LIBS})

find_package(GTest REQUIRED)
target_link_libraries(testColorConvert ${GTEST_BOTH_LIBRARIES})

find_package(RunningStats REQUIRED)
include_directories(${RUNNINGSTATS_INCLUDE_DIRS})

target_link_libraries(testColorConvert ${RUNNINGSTATS_LIBRARIES})
add_executable(plotcolorspace plotcolorspace.cpp)
target_link_libraries(plotcolorspace ${OpenCV_LIBS})
target_link_libraries(plotcolorspace ColorConvert)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
add_executable(optimize-color-circle optimize-color-circle.cpp)
target_link_libraries(optimize-color-circle
    ${CERES_LIBRARIES}
    ${OpenCV_LIBS}
    ColorConvert
    )

add_executable(testcolorflow testcolorflow.cpp)
target_link_libraries(testcolorflow
    ColorConvert
    ${GTEST_BOTH_LIBRARIES}
    )

add_executable(color-flow color-flow.cpp)
target_link_libraries(color-flow
    ColorConvert
    ${RUNNINGSTATS_LIBRARIES}
    )
