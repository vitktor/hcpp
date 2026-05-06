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
        Polynomial(
            std::vector<T> coefficients,
            Eigen::MatrixXi exponents,
            std::vector<Variable> variables
        )
            : coeffs(std::move(coefficients)), exps(std::move(exponents)), vars(std::move(variables))
            {
                
            }

        const std::vector<T>& getCoefficients() const
        {
            return coeffs;
        }

        void setCoefficients(std::vector<T> coefficients)
        {
            coeffs = std::move(coefficients);
        }

        const Eigen::MatrixXi& getExponents() const
        {
            return exps;
        }

        void setExponents(Eigen::MatrixXi exponents)
        {
            exps = std::move(exponents);
        }

        const std::vector<Variable>& getVariables() const
        {
            return vars;
        }

        void setVariables(std::vector<Variable> variables)
        {
            vars = std::move(variables);
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