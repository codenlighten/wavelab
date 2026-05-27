cmake_minimum_required(VERSION 3.24)
project(wavelab
    VERSION 0.0.1
    DESCRIPTION "Wave-based molecular geometry engine (FDTD)"
    LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Place all build artifacts under build/ in predictable locations.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ------------------------------------------------------------------
# Options
# ------------------------------------------------------------------
option(WAVELAB_BUILD_TESTS       "Build unit + integration tests"          ON)
option(WAVELAB_BUILD_GUI         "Build wavelab GUI (needs X11/GL dev libs)" ON)
option(WAVELAB_BUILD_CLI         "Build headless wavecli runner"            ON)
option(WAVELAB_BUILD_PYTHON      "Build pybind11 Python module"             ON)
option(WAVELAB_BUILD_BENCH       "Build perf benchmarks"                    OFF)
option(WAVELAB_DOUBLE_PRECISION  "Use double instead of float for fields"   OFF)
option(WAVELAB_USE_OPENMP        "Enable OpenMP parallelism"                ON)
option(WAVELAB_SANITIZE          "Enable ASan + UBSan in Debug builds"      OFF)
option(WAVELAB_WERROR            "Treat warnings as errors"                 ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(CompilerFlags)
include(Dependencies)

# ------------------------------------------------------------------
# Core library
# ------------------------------------------------------------------
add_subdirectory(src)

# ------------------------------------------------------------------
# Apps
# ------------------------------------------------------------------
add_subdirectory(apps)

# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------
if(WAVELAB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(WAVELAB_BUILD_PYTHON)
    add_subdirectory(python)
endif()

if(WAVELAB_BUILD_BENCH)
    add_subdirectory(bench)
endif()

# ------------------------------------------------------------------
# Build summary
# ------------------------------------------------------------------
message(STATUS "")
message(STATUS "wavelab ${PROJECT_VERSION} configured:")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler          : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Precision         : ${WAVELAB_PRECISION_NAME}")
message(STATUS "  OpenMP            : ${WAVELAB_OPENMP_STATUS}")
message(STATUS "  Tests             : ${WAVELAB_BUILD_TESTS}")
message(STATUS "  GUI app           : ${WAVELAB_BUILD_GUI}")
message(STATUS "  CLI app           : ${WAVELAB_BUILD_CLI}")
message(STATUS "  Python module     : ${WAVELAB_BUILD_PYTHON}")
message(STATUS "  Benchmarks        : ${WAVELAB_BUILD_BENCH}")
message(STATUS "  Sanitizers        : ${WAVELAB_SANITIZE}")
message(STATUS "")
