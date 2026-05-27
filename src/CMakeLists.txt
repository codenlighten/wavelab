# wavelab core library — interface-style API plus the .cpp implementations
# (PDB / SDF parsers today; FDTD instantiations could land here later).

add_library(wavelab_core STATIC
    molecule/parser_pdb.cpp
    molecule/parser_sdf.cpp
    score/spectral.cpp
    score/calibration.cpp
    io/fingerprint.cpp)
add_library(wavelab::core ALIAS wavelab_core)

target_include_directories(wavelab_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(wavelab_core PUBLIC
    wavelab::flags
    spdlog::spdlog
    Eigen3::Eigen
    pocketfft
    tomlplusplus::tomlplusplus)

# CUDA backend — separate library so wavelab_core stays nvcc-free.
if(WAVELAB_BUILD_CUDA)
    add_library(wavelab_cuda STATIC
        fdtd/fdtd_cuda.cu)
    add_library(wavelab::cuda ALIAS wavelab_cuda)

    target_include_directories(wavelab_cuda PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

    target_link_libraries(wavelab_cuda PUBLIC
        wavelab::core
        CUDA::cudart)

    # nvcc trips many of our strict host-compiler warnings on its own
    # generated code; keep the CUDA TU on relaxed warning settings.
    set_target_properties(wavelab_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
endif()
