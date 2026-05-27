# Centralized compiler flags. Apply by linking against wavelab::flags.

add_library(wavelab_flags INTERFACE)
add_library(wavelab::flags ALIAS wavelab_flags)

if(MSVC)
    target_compile_options(wavelab_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/W4 /permissive- /Zc:__cplusplus /Zc:preprocessor /utf-8>)
    if(WAVELAB_WERROR)
        target_compile_options(wavelab_flags INTERFACE
            $<$<COMPILE_LANGUAGE:CXX>:/WX>)
    endif()
else()
    # Scope to CXX so nvcc (which parses its own subset of Werror values)
    # doesn't get host-compiler flags shoved down its throat.
    target_compile_options(wavelab_flags INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic
        -Wshadow -Wconversion -Wsign-conversion
        -Wnon-virtual-dtor -Woverloaded-virtual
        -Wnull-dereference -Wdouble-promotion
        -Wformat=2>)
    if(WAVELAB_WERROR)
        target_compile_options(wavelab_flags INTERFACE
            $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
    endif()
endif()

# Precision selector
if(WAVELAB_DOUBLE_PRECISION)
    target_compile_definitions(wavelab_flags INTERFACE WAVELAB_DOUBLE_PRECISION=1)
    set(WAVELAB_PRECISION_NAME "double")
else()
    target_compile_definitions(wavelab_flags INTERFACE WAVELAB_DOUBLE_PRECISION=0)
    set(WAVELAB_PRECISION_NAME "float")
endif()

# OpenMP
if(WAVELAB_USE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(wavelab_flags INTERFACE OpenMP::OpenMP_CXX)
        target_compile_definitions(wavelab_flags INTERFACE WAVELAB_HAVE_OPENMP=1)
        set(WAVELAB_OPENMP_STATUS "enabled (${OpenMP_CXX_VERSION})")
    else()
        message(WARNING "OpenMP requested but not found; building serial.")
        set(WAVELAB_OPENMP_STATUS "requested but missing")
    endif()
else()
    set(WAVELAB_OPENMP_STATUS "disabled")
endif()

# Sanitizers (Debug only)
if(WAVELAB_SANITIZE AND NOT MSVC)
    target_compile_options(wavelab_flags INTERFACE
        $<$<CONFIG:Debug>:-fsanitize=address,undefined -fno-omit-frame-pointer>)
    target_link_options(wavelab_flags INTERFACE
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>)
endif()
