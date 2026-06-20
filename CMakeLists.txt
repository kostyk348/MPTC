cmake_minimum_required(VERSION 3.16)
project(mptc
    VERSION 1.0.0
    DESCRIPTION "Multi-Phase Ternary Computing library — branchless fixed-point math for low-power MCUs"
    LANGUAGES C)

# --------------------------------------------------------------------
#  Options
# --------------------------------------------------------------------
option(MPTC_BUILD_TESTS    "Build unit tests"              ON)
option(MPTC_BUILD_BENCH    "Build PC benchmark"            ON)
option(MPTC_USE_FLOAT      "Expose float32 overloads"      OFF)
option(MPTC_STRICT_BRANCHLESS "Forbid branches in hot path" ON)
option(MPTC_BUILD_DOCS     "Build Doxygen documentation"   OFF)

# --------------------------------------------------------------------
#  Library target
# --------------------------------------------------------------------
add_library(mptc STATIC
    src/mptc.c
    src/sdft.c
    src/integral_img.c)
add_library(mptc::mptc ALIAS mptc)

target_include_directories(mptc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_compile_features(mptc PUBLIC c_std_99)

# Math library link (for sqrt, cos, sin used in twiddle table generation).
target_link_libraries(mptc PUBLIC m)

if(MPTC_USE_FLOAT)
    target_compile_definitions(mptc PUBLIC MPTC_USE_FLOAT=1)
endif()

if(MPTC_STRICT_BRANCHLESS)
    target_compile_definitions(mptc PUBLIC MPTC_BRANCHLESS_STRICT=1)
endif()

# Recommended flags: aggressive optimization, no fast-math (we need IEEE here).
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(mptc PRIVATE
        -Wall -Wextra -Wno-unused-parameter -Wno-unused-function
        -O2 -fno-math-errno)
endif()

# --------------------------------------------------------------------
#  Tests
# --------------------------------------------------------------------
if(MPTC_BUILD_TESTS)
    enable_testing()
    add_executable(test_basic tests/test_basic.c)
    target_link_libraries(test_basic PRIVATE mptc)
    target_compile_definitions(test_basic PRIVATE MPTC_TEST_HOST=1)
    add_test(NAME test_basic COMMAND test_basic)
endif()

# --------------------------------------------------------------------
#  Benchmark
# --------------------------------------------------------------------
if(MPTC_BUILD_BENCH)
    add_executable(pc_benchmark examples/pc/pc_benchmark.c)
    target_link_libraries(pc_benchmark PRIVATE mptc)
    target_compile_definitions(pc_benchmark PRIVATE MPTC_TEST_HOST=1)
endif()

# --------------------------------------------------------------------
#  Install
# --------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS mptc
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# --------------------------------------------------------------------
#  Doxygen
# --------------------------------------------------------------------
if(MPTC_BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        add_custom_target(docs
            COMMAND ${CMAKE_COMMAND} -E copy ${DOXYGEN_IN} ${DOXYGEN_OUT}
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating Doxygen docs"
            VERBATIM)
    endif()
endif()
