#pragma once
#include "utils.hpp"
#include "chapter02.hpp"

struct ArnoldResult{
    Matrix* v;
    Matrix H;
    ArnoldResult(Matrix* v, Matrix H): v(v), H(H){};
    ~ArnoldResult(){
        delete[] v;
    };
};

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

    // 阿诺尔迪过程
    ArnoldResult ArnoldiProcess(const Matrix& r, int m);
    // 阿诺尔迪迭代法(循环型)
    Matrix ArnoldiSolve(int m = 10);

    // 广义极小残余算法(循环型)
    Matrix GMRESSolve(int m = 10);

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
    static void Demo3(int m);
    static void Demo31(int m);
};