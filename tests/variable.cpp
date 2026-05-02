#include <hc/core/variable.hpp>
#include <cassert>
#include <iostream>

int main() {
    auto vars = hc::make_variables({"x", "y", "z"});
    assert(vars.size() == 3);
    assert(vars[0].name == "x");
    assert(vars[1].name == "y");
    assert(vars[2].name == "z");

    hc::Variable a("a"), b("b");
    assert(a.name == "a");
    assert(b.name == "b");

    hc::Variable alpha("α"), beta("β");
    assert(alpha.name == "α");
    assert(beta.name == "β");

    std::cout << "All tests passed\n";
    return 0;
}
