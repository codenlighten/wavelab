add_executable(wavelab
    main.cpp
    ${CMAKE_SOURCE_DIR}/src/viz/view_1d.cpp
    ${CMAKE_SOURCE_DIR}/src/viz/view_2d.cpp)

target_link_libraries(wavelab PRIVATE
    wavelab::core
    wavelab::flags
    raylib
    imgui
    rlimgui)
