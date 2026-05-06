#include <vector>

#include <hc/core/polynomial.hpp>

namespace hc
{
    Polynomial<double> operator+(const Variable& lv, const Variable& rv)
    {
        if (lv == rv)
        {
            Eigen::MatrixXi exps = Eigen::MatrixXi::Identity(1, 1);
            return Polynomial<double>({2.0}, exps, {lv});
        }
        else
        {
            std::vector<Variable> vars;
            if (lv < rv)
                vars = {lv, rv};
            else
                vars = {rv, lv};
            
            Eigen::MatrixXi exps = Eigen::MatrixXi::Identity(2, 2);
            return Polynomial<double>({1.0, 1.0}, exps, vars);
        }
    }
}