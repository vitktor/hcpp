#pragma once
#include <vector>
#include <Eigen/Dense>

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

        void setExponents(Eigen:MatrixXi exponents)
        {
            exps = std::move(exponents);
        }

    private:
        std::vector<T> coeffs;
        Eigen::MatrixXi exps;
        std::vector<Variable> vars;
    };
}