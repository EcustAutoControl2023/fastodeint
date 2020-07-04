cmake_minimum_required(VERSION 3.8)
project(odeint_11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

add_subdirectory(pybind11)
pybind11_add_module(odeint_11 src/main.cpp)

target_link_libraries(odeint_11 PUBLIC OpenMP::OpenMP_CXX PUBLIC Boost::headers)
