add_executable(univariate_quadratic univariate_quadratic/main.cpp)
target_link_libraries(univariate_quadratic PRIVATE hc)
target_compile_definitions(univariate_quadratic PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/univariate_quadratic")

add_executable(circle_and_line circle_and_line/main.cpp)
target_link_libraries(circle_and_line PRIVATE hc)
target_compile_definitions(circle_and_line PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/circle_and_line")
