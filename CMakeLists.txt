cmake_minimum_required(VERSION 3.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(Reconstruction LANGUAGES CXX)

add_subdirectory(utils)

file(GLOB SRC_FILES *.cpp)

add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${CMAKE_PROJECT_NAME} Utils)
target_compile_definitions(${CMAKE_PROJECT_NAME} PUBLIC _USE_MATH_DEFINES)
