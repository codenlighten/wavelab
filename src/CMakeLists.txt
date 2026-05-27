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
