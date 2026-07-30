cmake_minimum_required(VERSION 3.10.2)
project(REF_version_werner)

set(CMAKE_BUILD_TYPE Release)

#Bring the headers, such as Student.h into the project
include_directories(include)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp" "src/*.c" "main.cpp")

add_executable(REF_version_werner.exe ${SOURCES})
