include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SCIP_INCLUDE_DIRS}
  ${SoPlex_INCLUDE_DIRS}
)

link_directories(
  ${SCIP_LIBRARY_DIRS}
  ${CMAKE_BINARY_DIR}/cpm/
)


add_executable(nonnegative-rank-bounds
  main.cpp
  scip_oracle.cpp
  enum_oracle.cpp
  slackmatrix.cpp
)

add_dependencies(nonnegative-rank-bounds cpm)

target_link_libraries(nonnegative-rank-bounds
  ${SCIP_LIBRARIES}
  cpm
  -lm
)
