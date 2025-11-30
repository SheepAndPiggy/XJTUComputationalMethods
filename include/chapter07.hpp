#pragma once
#include "utils.hpp"
#include <functional>

class NonlinearEquationSolver{
    public:
    NonlinearEquationSolver(std::function<double(double)> func);

    // 二分法
    double bisectionSolver(double min_x, double max_x);

    // 简单迭代法
    double simpleInterSolver(std::function<double(double)> phi, double x0 = 0);

    // 牛顿迭代法
    double newtonInterSolver(std::function<double(double)> dfunc, double x0 = 0);

    // 弦割法
    double secantInterSolver(double min_x, double max_x);

    // 斯特芬森加速技术
    double stephensonInterSolver(std::function<double(double)> phi, double x0 = 0);

    // 利用二分法找到解存在区间
    std::vector<std::pair<double, double>> findAllRootIntervals(
    double L, double R,
    double step);

    private:
    std::function<double(double)> func;
};

// 非线性方程组求解类：支持牛顿法、弦割法、布罗伊登法
class NonlinearSystemSolver {
    public:
    using Vector      = Matrix;                          // 约定：列向量用 n×1 的 Matrix 表示
    using Func        = std::function<Vector(const Vector&)>;     // F(x)
    using JacobianFun = std::function<Matrix(const Vector&)>;     // J(x)

    // 构造函数：传入非线性方程组 F(x)
    explicit NonlinearSystemSolver(Func func,
                                   double tol = 1e-8,
                                   int maxIter = 100);

    // 设置迭代精度和最大迭代次数
    void setTolerance(double t);
    void setMaxIter(int it);

    // 牛顿法：需要提供雅可比矩阵 J(x)
    Vector newtonSolver(const Vector& x0,
                        JacobianFun jacobian);

    // 弦割法（Chord Method）：使用固定雅可比矩阵（在初始点计算）
    Vector chordSolver(const Vector& x0,
                       JacobianFun jacobian);

    // 布罗伊登法（Broyden Method）：使用初始近似雅可比矩阵 B0
    Vector broydenSolver(const Vector& x0,
                         const Matrix& B0);

    private:
    Func   func;       // 非线性方程组 F(x)
    double tol;        // 收敛精度
    int    maxIter;    // 最大迭代次数

    // 求解线性方程组 A * dx = b （这里简单用 A^{-1} * b）
    Vector solveLinearSystem(const Matrix& A,
                             const Vector& b) const;

    // 计算列向量 a,b 的外积：result(i,j) = a(i,0) * b(j,0)
    Matrix outerProduct(const Vector& a,
                        const Vector& b) const;
};

struct Chapter07Demos{
    static void Demo1();
    static void Demo2();
};
