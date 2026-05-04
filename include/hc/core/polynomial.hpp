#pragma once
#include <vector>
#include <Eigen/Dense>

#include <hc/core/variable.hpp>

namespace hc
{
    template <typename T>
    class Polynomial
    {
    public:
        Polynomial(Eigen::MatrixXi exponents, std::vector<T> coefficients)
            : exps(std::move(exponents)), coeffs(std::move(coefficients)) {}

        const std::vector<T>& getCoefficients() const
        {
            return coeffs;
        }

        void setCoefficients(std::vector<T> coefficients)
        {
            coeffs = std::move(coefficients);
        }

        const std::vector<Variable>& getVariables() const
        {
            return vars;
        }

        void setVariables(std::vector<Variable> variables)
        {
            vars = std::move(variables);
        }

        const Eigen::MatrixXi& getExponents() const
        {
            return exps;
        }

        void setExponents(Eigen::MatrixXi exponents)
        {
            exps = std::move(exponents);
        }

        int degree() const
        {
            if (exps.cols() == 0) return -1;
            return static_cast<int>(exps.colwise().sum().maxCoeff());
        }

        int degree(const Variable& var) const
        {
            auto it = std::find(vars.begin(), vars.end(), var);
            if (it == vars.end()) return 0;
            auto id = it - vars.begin();
            return static_cast<int>(exps.row(id).maxCoeff());
        }

    private:
        std::vector<T> coeffs;
        Eigen::MatrixXi exps;
        std::vector<Variable> vars;
    };

Polynomial<double> operator+(const Variable& lv, const Variable& rv);
Polynomial<double> operator-(const Variable& lv, const Variable& rv);
Polynomial<double> operator*(const Variable& lv, const Variable& rv);
Polynomial<double> pow(const Variable& lv, int exp);

}