if(WAVELAB_BUILD_GUI)
    add_subdirectory(wavelab)
endif()

if(WAVELAB_BUILD_CLI)
    add_subdirectory(wavecli)
endif()
