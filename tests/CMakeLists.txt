include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)

# All unit tests share a single binary that owns doctest's main(), so we
# get one fast compile + one test invocation. Each test_*.cpp just adds
# TEST_CASEs into this binary.

set(WAVELAB_TEST_SOURCES
    unit/test_main.cpp
    unit/test_types.cpp
    unit/test_grid.cpp
    unit/test_field.cpp
    unit/test_stepper_1d.cpp
    unit/test_energy_1d.cpp
    unit/test_similarity.cpp
    unit/test_fdtd_2d.cpp
    unit/test_entropy.cpp
    unit/test_spectral.cpp
    unit/test_fingerprint.cpp
    unit/test_fdtd_3d.cpp
    integration/test_1ubq_pipeline.cpp
    integration/test_1ubq_full_3d.cpp
    integration/test_determinism.cpp
    integration/test_binder_discrimination.cpp)

if(WAVELAB_BUILD_CUDA)
    list(APPEND WAVELAB_TEST_SOURCES unit/test_fdtd_cuda.cpp)
endif()

add_executable(wavelab_tests ${WAVELAB_TEST_SOURCES})

target_link_libraries(wavelab_tests PRIVATE
    wavelab::core
    wavelab::flags
    doctest::doctest)

if(WAVELAB_BUILD_CUDA)
    target_link_libraries(wavelab_tests PRIVATE wavelab::cuda)
endif()

# Auto-register every TEST_CASE with CTest.
doctest_discover_tests(wavelab_tests)

# Python smoke test (only when the binding module is built).
if(WAVELAB_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_FOUND)
        add_test(NAME wavelab_python_smoke
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/python/test_pybind.py)
        set_tests_properties(wavelab_python_smoke PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endif()
endif()
