#include <hc/core/variable.hpp>

namespace hc
{
    std::vector<Variable> make_variables(std::initializer_list<std::string> names)
    {
        std::vector<Variable> vars;
        for (const auto& name : names)
            vars.emplace_back(name);
        return vars;
    }
}