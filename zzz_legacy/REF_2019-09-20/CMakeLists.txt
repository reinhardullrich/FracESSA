cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type to Release
set(CMAKE_BUILD_TYPE Release)

project(REF)

include_directories(include)

if(UNIX AND NOT APPLE)
    include_directories(/home/reinhard/REF/boost_1_71_0_linux)
    add_executable(REF_linux64.exe src/argtable3.c src/EssFinder.cpp src/main.cpp)
    # Set maximum optimization flags with all optimizations on, -pg off
    target_compile_options(REF_linux64.exe PRIVATE 
        -O3 
        -march=native 
        -flto 
        -DNDEBUG
    )
    target_link_options(REF_linux64.exe PRIVATE -flto)
endif()

if(WIN32 OR MSVC OR MSYS OR MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
    include_directories(/home/reinhard/REF/boost_1_71_0_win)
    add_executable(REF_win64 src/argtable3.c src/EssFinder.cpp src/main.cpp)
endif()
