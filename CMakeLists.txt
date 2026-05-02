cmake_minimum_required(VERSION 3.20)
project(hcpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

add_library(hc
    src/variable.cpp
    src/polynomial.cpp
    src/system.cpp
    src/differentiation.cpp
    src/tracker.cpp
    src/corrector.cpp
    src/runge_kutta.cpp
    src/step_control.cpp
    src/straight_line_homotopy.cpp
    src/parameter_homotopy.cpp
    src/total_degree.cpp
    src/solver.cpp
    src/cauchy_endgame.cpp
)

target_include_directories(hc PUBLIC include)
target_link_libraries(hc PUBLIC Eigen3::Eigen)

add_executable(tests tests/variable.cpp)
target_link_libraries(tests PRIVATE hc)
