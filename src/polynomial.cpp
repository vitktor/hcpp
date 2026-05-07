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

    Polynomial<double> operator-(const Variable& lv, const Variable& rv)
    {
        if (lv == rv)
            return Polynomial<double>(0.0);
        else
        {
            Eigen::MatrixXi exps = Eigen::MatrixXi::Identity(2, 2);
            if (lv < rv)
                return Polynomial<double>({1.0, -1.0}, exps, {lv, rv});
            else
                return Polynomial<double>({-1.0, 1.0}, exps, {rv, lv});
        }
    }

    Polynomial<double> operator*(const Variable& lv, const Variable& rv)
    {
        if (lv == rv)
        {
            Eigen::MatrixXi exps(1,1);
            exps(0,0) = 2;
            return Polynomial<double>({1.0}, exps, {lv});
        }
        else
        {
            std::vector<Variable> vars;
            if (lv < rv)
                vars = {lv, rv};
            else
                vars = {rv, lv};

            Eigen::MatrixXi exps(2,1);
            exps << 1, 1;
            return Polynomial<double>({1.0}, exps, vars);
        }
    }

    Polynomial<double> pow(const Variable& var, int exp)
    {
        if (exp == 0)
            return Polynomial<double>(1.0);
        Eigen::MatrixXi exps(1,1);
        exps(0,0) = exp;
        return Polynomial<double>({1.0}, exps, {var});
    }
}