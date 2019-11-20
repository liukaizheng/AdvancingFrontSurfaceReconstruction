file(GLOB SRC_FILES *.cpp)

add_library(GlpkSolver OBJECT ${SRC_FILES})

if(NOT GLPK_ROOT)
    set(GLPK_ROOT "D:/ProgramFiles/glpk-4.65")
endif()

target_include_directories(GlpkSolver PRIVATE ${GLPK_ROOT}/src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_library(GLPK_LIB NAMES glpk PATHS ${GLPK_ROOT}/w64)
if(NOT GLPK_LIB)
    message(FATAL_ERROR "cananot find glpk library in ${GLPK_LIB}/w64")
endif()

target_link_libraries(GlpkSolver PUBLIC ${GLPK_LIB})
