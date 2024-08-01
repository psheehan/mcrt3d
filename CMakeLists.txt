cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(Kokkos REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Python COMPONENTS Interpreter NumPy REQUIRED)

pybind11_add_module(cpu src/pinballrt.cc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cpu PUBLIC Kokkos::kokkos)
target_link_libraries(cpu PUBLIC Python::NumPy)

install(TARGETS cpu LIBRARY DESTINATION pinballrt)
