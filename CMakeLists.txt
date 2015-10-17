cmake_minimum_required(VERSION 2.8)
project(kmc)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Executables:
add_executable(a main.cpp bc.cpp boffsetassign.cpp dipassign.cpp getenergy.cpp
realsum.cpp recsum.cpp selfint.cpp structure2.cpp)
