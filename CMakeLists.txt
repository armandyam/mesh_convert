cmake_minimum_required(VERSION 3.0.0)

project(myproject LANGUAGES CXX)

set(Omega_h_PREFIX "" CACHE STRING "Location of Omega_h install")

find_package(Omega_h 3.0.1 REQUIRED PATHS "${Omega_h_PREFIX}")

add_executable(omega_h netgen2omega.cpp)
target_link_libraries(omega_h Omega_h::omega_h)
