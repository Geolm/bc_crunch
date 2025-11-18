cmake_minimum_required(VERSION 3.25)
project(bc_crunch)

if(CMAKE_GENERATOR STREQUAL "Xcode")
    set(CMAKE_MACOSX_BUNDLE TRUE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_C_STANDARD 11)
if(MSVC)
    set(CMAKE_C_FLAGS_DEBUG "/W3 /Zi")
else()
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG -Wall -Wextra ")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Wextra")
endif()

add_executable(test
    ./test/test.c
    ./bc_crunch.c
)

if(UNIX AND NOT APPLE)
    target_link_libraries(test m)
endif()


add_executable(benchmark
    ./test/benchmark.c
    ./bc_crunch.c
)

if(UNIX AND NOT APPLE)
    target_link_libraries(benchmark m)
endif()