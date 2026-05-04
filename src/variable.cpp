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

    MergedVars merge_vars(
        const std::vector<Variable>& vars1,
        const std::vector<Variable>& vars2
    )
    {
        MergedVars result;
        result.vars.reserve(vars1.size() + vars2.size());
        result.map1.resize(vars1.size());
        result.map2.resize(vars2.size());

        size_t i = 0, j = 0, k = 0;
        while (i < vars1.size() && j < vars2.size())
        {
            if (vars1[i] < vars2[j])
            {
                result.map1[i++] = k;
                result.vars.push_back(vars1[i-1]);
            }
            else if (vars2[j] < vars1[i])
            {
                result.map2[j++] = k;
                result.vars.push_back(vars2[j-1]);
            }
            else
            {
                result.map1[i++] = result.map2[j++] = k;
                result.vars.push_back(vars1[i-1]);
            }
            ++k;
        }
        while (i < vars1.size()) { result.map1[i++] = k++; result.vars.push_back(vars1[i-1]); }
        while (j < vars2.size()) { result.map2[j++] = k++; result.vars.push_back(vars2[j-1]); }

        return result;
    }
}