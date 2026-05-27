add_executable(bench_fdtd bench_fdtd.cpp)

target_link_libraries(bench_fdtd PRIVATE
    wavelab::core
    wavelab::flags)

if(WAVELAB_BUILD_CUDA)
    target_link_libraries(bench_fdtd PRIVATE wavelab::cuda)
    target_compile_definitions(bench_fdtd PRIVATE WAVELAB_BENCH_HAVE_CUDA=1)
else()
    target_compile_definitions(bench_fdtd PRIVATE WAVELAB_BENCH_HAVE_CUDA=0)
endif()
