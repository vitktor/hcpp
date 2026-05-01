#include <hc/core/variable.hpp>
#include <iostream>

int main() {
    auto [x, y, z] = hc::make_variables("x", "y", "z");
    std::cout << x.name << "\n";
    std::cout << y.name << "\n";
    std::cout << z.name << "\n";

    auto vars = hc::make_variables({"a", "b", "c"});
    for (const auto& v : vars)
        std::cout << v.name << "\n";

    return 0;
}
