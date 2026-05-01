#pragma once
#include <string>
#include <tuple>
#include <vector>

namespace hc
{
    struct Variable
    {
        std::string name;
        Variable(std::string name) : name(std::move(name)) {}
    };

    template <typename... Names>
    auto make_variables(Names&&... names)
    {
        return std::make_tuple(Variable(std::forward<Names>(names))...);
    }

    inline std::vector<Variable> make_variables(std::initializer_list<std::string> names)
    {
        std::vector<Variable> vars;
        for (const auto& name : names)
            vars.emplace_back(name);
        return vars;
    }
}
