add_executable(tracking_univariate_quadratic tracking/univariate_quadratic/main.cpp)
target_link_libraries(tracking_univariate_quadratic PRIVATE hc)
target_compile_definitions(tracking_univariate_quadratic PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tracking/univariate_quadratic")

add_executable(tracking_circle_and_line tracking/circle_and_line/main.cpp)
target_link_libraries(tracking_circle_and_line PRIVATE hc)
target_compile_definitions(tracking_circle_and_line PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tracking/circle_and_line")

add_executable(solving_univariate_quadratic solving/univariate_quadratic/main.cpp)
target_link_libraries(solving_univariate_quadratic PRIVATE hc)
target_compile_definitions(solving_univariate_quadratic PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/solving/univariate_quadratic")

add_executable(solving_circle_and_parabola solving/circle_and_parabola/main.cpp)
target_link_libraries(solving_circle_and_parabola PRIVATE hc)
target_compile_definitions(solving_circle_and_parabola PRIVATE
    HC_EXAMPLE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/solving/circle_and_parabola")
