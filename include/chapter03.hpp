#pragma once
#include "utils.hpp"

class IterationSolver{
    public:
    Matrix D, E, F, D_inv;  // 声明DEF矩阵和D的逆矩阵

    // 构造函数，生成D、E、F
    IterationSolver(const Matrix& A, const Matrix& b, unsigned int max_iter = 10000, double epsilon = 1e-8);

    // 雅可比迭代法
    Matrix jacobiSolve();
    Matrix SORSolve(double w);
    Matrix gaussSolve();

    // 共轭梯度法
    Matrix ConjugateGradientSolve();

    // 误差函数
    double convergenceError(const Matrix& x_old, const Matrix& x_new);
    double residualNormError(const Matrix& x);

    private:
    Matrix A;
    Matrix b;
    unsigned int max_iter;
    double epsilon;

    // 基础迭代法
    Matrix baseSolve(const Matrix& B, const Matrix& g);
};

// 测试用例
struct Chapter03Demos{
    static void Demo1();
    static void Demo2();
};