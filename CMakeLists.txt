cmake_minimum_required(VERSION 3.10)
project(math)

set(CMAKE_CXX_STANDARD 17)

include_directories("include")
file(GLOB_RECURSE SRC_FILE ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILE})