cmake_minimum_required(VERSION 3.0.0)

project(myproject LANGUAGES CXX)

set(Omega_h_PREFIX "" CACHE STRING "Location of Omega_h install")

find_package(Omega_h 3.2.0 REQUIRED PATHS "${Omega_h_PREFIX}")

add_executable(mesh_read netgen_read.cpp)
target_link_libraries(mesh_read Omega_h::omega_h)
