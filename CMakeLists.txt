cmake_minimum_required(VERSION 3.0)

project("Nonnegative rank bounds" VERSION 0.1)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_definitions("-std=gnu++11")

find_package(SCIP REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(cpm)
add_subdirectory(src)


