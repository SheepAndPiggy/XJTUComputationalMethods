#pragma once
#include "utils.hpp"
#include "chapter03.hpp"
#include <functional>

class LeastSquaresApproximator{
    public:
    LeastSquaresApproximator(Matrix x, Matrix y);

    // 构造形如gk=x^k基底的最小二乘拟合函数
    std::function<double(double)> standardPolynomial(unsigned int n);
    // 构造勒让德多项式组成的正交基底的最小二乘函数
    std::function<double(double)> legendrePolynomial(unsigned int n);

    private:
    Matrix x;
    Matrix y;
    int n;
    double max_x;
    double min_x;

    // 内积计算函数，计算任意两个函数g1,g2在点集上的w内积
    double innerProduct(std::function<double(double)> g1, std::function<double(double)> g2,
    std::function<double(double)> w = nullptr);

    // 基本多项式基函数生成函数
    std::function<double(double)> standardPolynomialBase(int k);
    // 由勒让德多项式的三项递推关系得出的基函数生成函数
    std::function<double(double)> legendrePolynomialBase(int k, std::function<double(double)> g1, std::function<double(double)> g2);

    // 基函数相互正交的拟合函数：适用于一组正交基底
    std::function<double(double)> orthogonalBasisFit(std::vector<std::function<double(double)>> funcs, 
    std::function<double(double)> w = nullptr);

    // 通过求解正规方程组的拟合函数：适用于一组线性无关但不正交的基底
    std::function<double(double)> normalEquationsFit(std::vector<std::function<double(double)>> funcs, 
    std::function<double(double)> w = nullptr);
};

// 测试用例
struct Chapter05Demos{
    static void Demo1();
};