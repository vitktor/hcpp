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
        bool operator==(const Variable& other) const { return name == other.name; }
        bool operator<(const Variable& other) const { return name < other.name; }
    };

    struct MergedVars
    {
        std::vector<Variable> vars;
        std::vector<int> map1;
        std::vector<int> map2;
    };

    std::vector<Variable> make_variables(std::initializer_list<std::string> names);
    MergedVars merge_sorted_vars(
        const std::vector<Variable>& vars1,
        const std::vector<Variable>& vars2
    );
}
