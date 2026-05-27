# Python binding module: imported as `import wavelab` once the build
# directory is on PYTHONPATH.

pybind11_add_module(wavelab MODULE wavelab_py.cpp)

target_link_libraries(wavelab PRIVATE
    wavelab::core
    wavelab::flags)

# Drop the lib prefix so the file lands as `wavelab.<ext>` regardless of
# platform — pybind11 normally handles this, but we override OUTPUT_NAME
# explicitly in case the CMake module dir was set elsewhere.
set_target_properties(wavelab PROPERTIES
    OUTPUT_NAME wavelab
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
