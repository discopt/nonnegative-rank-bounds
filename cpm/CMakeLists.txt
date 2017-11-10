include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SoPlex_INCLUDE_DIRS}
)

link_directories(
  ${SoPlex_LIBRARY_DIRS}
)

add_library(cpm
  core.cpp
  separation_oracle.cpp
  solver.cpp
  solver_soplex.cpp
)

target_link_libraries(cpm
  ${SoPlex_LIBRARIES}
  -lm
)

install(
  TARGETS cpm
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  COMPONENT library
)

install(
  DIRECTORY .
  DESTINATION include/cpm
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h"
)
