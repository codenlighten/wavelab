add_executable(wavecli main.cpp)

target_link_libraries(wavecli PRIVATE
    wavelab::core
    wavelab::flags)
