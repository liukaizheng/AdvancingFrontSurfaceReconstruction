file(GLOB UTILS_FILE *.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../eigen3/cmake)
message("CMKAE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
find_package(Eigen3 REQUIRED)

add_subdirectory(tetgen)
add_library(Utils OBJECT ${UTILS_FILE})

target_link_libraries(Utils PRIVATE tet)
target_include_directories(Utils PUBLIC ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
