#include "chapter07.hpp"
#include <functional>

NonlinearEquationSolver::NonlinearEquationSolver(std::function<double(double)> func): func(func){};

std::vector<std::pair<double, double>> NonlinearEquationSolver::findAllRootIntervals(
    double L, double R,
    double step)
{
    std::vector<std::pair<double, double>> intervals;
    auto f = this->func;
    double x_left = L;
    double f_left = f(x_left);

    for (double x_right = L + step; x_right <= R; x_right += step) {

        double f_right = f(x_right);

        // 检查是否跨过 x 轴（异号）
        if (f_left * f_right < 0) {
            intervals.push_back({x_left, x_right});
        }

        // 移动窗口
        x_left = x_right;
        f_left = f_right;
    }

    return intervals;
}

double NonlinearEquationSolver::bisectionSolver(double min_x, double max_x){
    double epsilon1, epsilon2;
    epsilon1 = epsilon2 = 1e-8;

    double a = min_x;
    double b = max_x;
    double x = (a + b) / 2;
    double sep = (b - a) / 2;
    double k = 1;
    while (std::abs(this->func(x)) > epsilon1 && sep > epsilon2){
        if (this->func(a) * this->func(x) < 0)
            b = x;
        else 
            a = x;
        x = (a + b) / 2;
        sep = (b - a) / 2;
        k += 1;
    }
    std::cout << "迭代法共用" << k << "步" << std::endl;
    return x;
}

double NonlinearEquationSolver::simpleInterSolver(std::function<double(double)> phi, double x0){
    double epsilon = 1e-8;
    double max_num = 1e4;

    double x = x0;
    double k = 1;
    while (std::abs(this->func(x)) > epsilon && k < max_num){
        x = phi(x);
        k += 1;
    }
    std::cout << "迭代法共用" << k << "步" << std::endl;
    return x;
}

double NonlinearEquationSolver::newtonInterSolver(std::function<double(double)> dfunc, double x0){
    auto phi = [this, dfunc](double x){
        return x - this->func(x) / dfunc(x);
    };
    return this->simpleInterSolver(phi, x0);
}

double NonlinearEquationSolver::secantInterSolver(double min_x, double max_x){
    double x1 = min_x;
    double x2 = max_x;
    double epsilon = 1e-8;

    int k = 1;
    while (std::abs(this->func(x2)) > epsilon){
        double y1 = this->func(x1);
        double y2 = this->func(x2);
        double temp = x2;
        x2 = x1 - y1 * (x2 - x1) / (y2 - y1);
        x1 = temp;
        k += 1;
    }
    std::cout << "迭代法共用" << k << "步" << std::endl;
    return x2;
}

double NonlinearEquationSolver::stephensonInterSolver(std::function<double(double)> phi, double x0){
    double x1 = phi(x0);
    double x2 = phi(x1);
    double x_ = x2 - std::pow(x2 - x1, 2) / (x2 + x0 - 2 * x1);
    double epsilon = 1e-8;

    int k = 1;
    while (std::abs(x_ - x2) > epsilon){
        x1 = phi(x_);
        x2 = phi(x1);
        x_ = x2 - std::pow(x2 - x1, 2) / (x2 + x_ - 2 * x1);
        k += 1;
    }
    std::cout << "迭代法共用" << k << "步" << std::endl;
    return x_;
}

NonlinearSystemSolver::NonlinearSystemSolver(Func func_,
                                             double tol_,
                                             int maxIter_)
    : func(std::move(func_)), tol(tol_), maxIter(maxIter_) {}

void NonlinearSystemSolver::setTolerance(double t) {
    tol = t;
}

void NonlinearSystemSolver::setMaxIter(int it) {
    maxIter = it;
}

NonlinearSystemSolver::Vector
NonlinearSystemSolver::solveLinearSystem(const Matrix& A,
                                         const Vector& b) const {
    // 简单做法：dx = A^{-1} * b
    // 如果你有更好的线性求解器（如高斯消元），可以替换这里
    Matrix Ainv = A.inv();
    return Ainv * b;
}

Matrix NonlinearSystemSolver::outerProduct(const Vector& a,
                                           const Vector& b) const {
    // 假设 a、b 都是 n×1 的列向量
    Matrix result(a.rows, b.rows);
    for (size_t i = 0; i < a.rows; ++i) {
        for (size_t j = 0; j < b.rows; ++j) {
            result(i, j) = a(i, 0) * b(j, 0);
        }
    }
    return result;
}

// ------------------------- 牛顿法 -------------------------

NonlinearSystemSolver::Vector
NonlinearSystemSolver::newtonSolver(const Vector& x0,
                                    JacobianFun jacobian) {
    Vector x = x0;

    for (int k = 0; k < maxIter; ++k) {
        Vector Fx = func(x);

        // 判断是否已经够精确
        if (norm(Fx) < tol) {
            // std::cout << "Newton converged in " << k << " iterations.\n";
            break;
        }

        Matrix J = jacobian(x);

        // 解 J * dx = -F(x)
        Vector dx = solveLinearSystem(J, -1.0 * Fx);

        x = x + dx;

        // 也可以用步长判断
        if (norm(dx) < tol) {
            break;
        }
    }

    return x;
}

// ------------------------- 弦割法（Chord Method） -------------------------

NonlinearSystemSolver::Vector
NonlinearSystemSolver::chordSolver(const Vector& x0,
                                   JacobianFun jacobian) {
    Vector x = x0;

    // 在初值处计算一次雅可比矩阵，并固定不再更新
    Matrix J0 = jacobian(x0);
    Matrix J0inv = J0.inv();

    for (int k = 0; k < maxIter; ++k) {
        Vector Fx = func(x);

        if (norm(Fx) < tol) {
            // std::cout << "Chord method converged in " << k << " iterations.\n";
            break;
        }

        // 用固定的 J0^{-1} 做更新：dx = - J0^{-1} F(x)
        Vector dx = J0inv * (-1.0 * Fx);

        x = x + dx;

        if (norm(dx) < tol) {
            break;
        }
    }

    return x;
}

// ------------------------- 布罗伊登法（Broyden Method） -------------------------

NonlinearSystemSolver::Vector
NonlinearSystemSolver::broydenSolver(const Vector& x0,
                                     const Matrix& B0) {
    Vector x = x0;
    Matrix B = B0;   // B 近似雅可比矩阵 J

    for (int k = 0; k < maxIter; ++k) {
        Vector Fx = func(x);

        if (norm(Fx) < tol) {
            // std::cout << "Broyden method converged in " << k << " iterations.\n";
            break;
        }

        // 解 B * dx = -F(x)
        Vector dx = solveLinearSystem(B, -1.0 * Fx);

        Vector x_new = x + dx;
        Vector Fx_new = func(x_new);

        Vector s = x_new - x;     // s_k = x_{k+1} - x_k
        Vector y = Fx_new - Fx;   // y_k = F_{k+1} - F_k

        // Broyden 更新：B_{k+1} = B_k + (y_k - B_k s_k) s_k^T / (s_k^T s_k)
        Vector Bs = B * s;        // B_k s_k
        Vector u  = y - Bs;       // y_k - B_k s_k

        double denom = norm(s);
        denom *= denom;           // denom = s^T s = ||s||^2

        if (denom > 0.0) {
            Matrix upd = outerProduct(u, s) / denom;
            B = B + upd;
        }

        x = x_new;
    }

    return x;
}


void Chapter07Demos::Demo1(){
    // 定义原函数、迭代格式、原函数导数的函数
    auto func = [](double x){
        return std::pow(x, 6) - 5 * std::pow(x, 5) + 3 * std::pow(x, 4) + 
        std::pow(x, 3) - 7 * std::pow(x, 2) + 7 * x - 20;
    };
    auto phi = [](double x){
        double val = 5 * std::pow(x, 5) - 3 * std::pow(x, 4)
                - std::pow(x, 3)     + 7 * std::pow(x, 2)
                - 7 * x + 20;
        return std::pow(val, 1.0 / 6.0);
    };
    auto dfunc = [](double x){
        return 6 * std::pow(x, 5)
            - 25 * std::pow(x, 4)
            + 12 * std::pow(x, 3)
            + 3 * std::pow(x, 2)
            - 14 * x
            + 7;
    };

    NonlinearEquationSolver solver(func);

    std::vector<std::pair<double, double>> ranges = solver.findAllRootIntervals(-1, 5, 0.01);
    double min_x = ranges[0].first;
    double max_x = ranges[0].second;

    std::cout << "简单";
    double x0 = solver.simpleInterSolver(phi, (min_x + max_x) / 2);
    std::cout << x0 << std::endl;

    std::cout << "牛顿";
    double x1 = solver.newtonInterSolver(dfunc, (min_x + max_x) / 2);
    std::cout << x1 << std::endl;

    std::cout << "弦割";
    double x2 = solver.secantInterSolver(min_x, max_x);
    std::cout << x2 << std::endl;

    std::cout << "使用斯蒂芬森加速技术的简单";
    double x3 = solver.stephensonInterSolver(phi, (min_x + max_x) / 2);
    std::cout << x3 << std::endl;
}

using Vector = Matrix;

void Chapter07Demos::Demo2() {

    /********************* (1) 三元非线性方程组 *************************
     * {
     *   x1^2 + x2^2 + x3^2 - 1.0 = 0
     *   2 x1^2 + x2^2 - 4 x3     = 0
     *   3 x1^2 - 4 x2^2 + x3^2   = 0
     * }
     *  给定初始向量 x^(0) = (1.0, 1.0, 1.0)^T
     *******************************************************************/
    auto F1 = [](const Vector& x) -> Vector {
        Vector fx(3, 1);
        double x1 = x(0, 0);
        double x2 = x(1, 0);
        double x3 = x(2, 0);

        fx(0, 0) = x1 * x1 + x2 * x2 + x3 * x3 - 1.0;
        fx(1, 0) = 2.0 * x1 * x1 + x2 * x2 - 4.0 * x3;
        fx(2, 0) = 3.0 * x1 * x1 - 4.0 * x2 * x2 + x3 * x3;
        return fx;
    };

    // 雅可比矩阵 J1(x)
    auto J1 = [](const Vector& x) -> Matrix {
        Matrix J(3, 3);
        double x1 = x(0, 0);
        double x2 = x(1, 0);
        double x3 = x(2, 0);

        // f1 = x1^2 + x2^2 + x3^2 - 1
        J(0, 0) = 2.0 * x1;
        J(0, 1) = 2.0 * x2;
        J(0, 2) = 2.0 * x3;

        // f2 = 2 x1^2 + x2^2 - 4 x3
        J(1, 0) = 4.0 * x1;
        J(1, 1) = 2.0 * x2;
        J(1, 2) = -4.0;

        // f3 = 3 x1^2 - 4 x2^2 + x3^2
        J(2, 0) = 6.0 * x1;
        J(2, 1) = -8.0 * x2;
        J(2, 2) = 2.0 * x3;

        return J;
    };

    // 初始向量
    Vector x0_1(3, 1);
    x0_1(0, 0) = 1.0;
    x0_1(1, 0) = 1.0;
    x0_1(2, 0) = 1.0;

    NonlinearSystemSolver solver1(F1, 1e-10, 50);

    Vector root1_newton   = solver1.newtonSolver(x0_1, J1);
    Vector root1_chord    = solver1.chordSolver(x0_1, J1);
    Matrix B0_1           = J1(x0_1);
    Vector root1_broyden  = solver1.broydenSolver(x0_1, B0_1);

    std::cout << "\n方程组 (1) 计算结果:\n";
    std::cout << "牛顿法:\n"   << root1_newton.transpose()  << "\n";
    std::cout << "弦割法:\n"   << root1_chord.transpose()  << "\n";
    std::cout << "布罗伊登法:\n"   << root1_broyden.transpose() << "\n";


    /********************* (2) 二元非线性方程组 *************************
     * {
     *   cos(x1^2 + 0.4 x2) + x1^2 + x2^2 - 1.6 = 0
     *   1.5 x1^2 - (1 / 0.36) x2^2 - 1.0       = 0
     * }
     *  给定初始向量 x^(0) = (1.04, 0.47)^T
     *******************************************************************/
    auto F2 = [](const Vector& x) -> Vector {
        Vector fx(2, 1);
        double x1 = x(0, 0);
        double x2 = x(1, 0);

        double g = x1 * x1 + 0.4 * x2;
        const double C = 1.0 / 0.36;   // 等于 2.777...

        fx(0, 0) = std::cos(g) + x1 * x1 + x2 * x2 - 1.6;
        fx(1, 0) = 1.5 * x1 * x1 - C * x2 * x2 - 1.0;

        return fx;
    };

    // 雅可比矩阵 J2(x)
    auto J2 = [](const Vector& x) -> Matrix {
        Matrix J(2, 2);
        double x1 = x(0, 0);
        double x2 = x(1, 0);

        double g = x1 * x1 + 0.4 * x2;
        double s = std::sin(g);
        const double C = 1.0 / 0.36;

        // f1_x1, f1_x2
        J(0, 0) = -2.0 * x1 * s + 2.0 * x1;  // 2x1(1 - sin(g))
        J(0, 1) = -0.4 * s + 2.0 * x2;

        // f2_x1, f2_x2
        J(1, 0) = 3.0 * x1;
        J(1, 1) = -2.0 * C * x2;

        return J;
    };

    Vector x0_2(2, 1);
    x0_2(0, 0) = 1.04;
    x0_2(1, 0) = 0.47;

    NonlinearSystemSolver solver2(F2, 1e-10, 50);

    Vector root2_newton  = solver2.newtonSolver(x0_2, J2);
    Vector root2_chord   = solver2.chordSolver(x0_2, J2);
    Matrix B0_2          = J2(x0_2);
    Vector root2_broyden = solver2.broydenSolver(x0_2, B0_2);

    std::cout << "\n方程组 (2) 计算结果:\n";
    std::cout << "牛顿法:\n"   << root2_newton.transpose()  << "\n";
    std::cout << "弦割法:\n"   << root2_chord.transpose()   << "\n";
    std::cout << "布罗伊登法:\n"   << root2_broyden.transpose() << "\n";
}
