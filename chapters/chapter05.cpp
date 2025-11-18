#include "utils.hpp"
#include "chapter05.hpp"
#include "chapter03.hpp"
#include <functional>

LeastSquaresApproximator::LeastSquaresApproximator(Matrix x, Matrix y):
x(x), y(y){
    if (x.rows != y.rows)
        throw std::invalid_argument("采样点位置和采样点取值维度必须相同！");
    this->n = x.rows;

    this->max_x = this->x(0, 0);
    this->min_x = this->x(0, 0);
    for (int i = 0; i < this->n; ++i){
        if (this->max_x < this->x(i, 0))
            this->max_x = this->x(i, 0);
        if (this->min_x > this->x(i, 0))
            this->min_x = this->x(i, 0);
    }
}

double LeastSquaresApproximator::innerProduct(std::function<double(double)> g1, std::function<double(double)> g2,
std::function<double(double)> w){
    double result = 0;
    // 如果没有传入w，则令所有w为1即不加权内积
    if (w == nullptr)
        w = [](double x){
            return 1;
        };

    for (int i = 0; i < this->n; ++i){
        double x = this->x(i, 0);
        if (g2 != nullptr)
            result += w(x) * g1(x) * g2(x);
        else
            result += w(x) * g1(x) * this->y(i, 0);
    }
    return result;
}

std::function<double(double)> LeastSquaresApproximator::standardPolynomialBase(int k){
    auto func = [k](double x){
        return std::pow(x, k);
    };
    return func;
}

std::function<double(double)> LeastSquaresApproximator::legendrePolynomialBase(int k, std::function<double(double)> g1, std::function<double(double)> g2){
    auto func = [this, k, g1, g2](double x){
        double t = 2 * (x - this->min_x) / (this->max_x - this->min_x) - 1;
        return (2 * k - 1) / (double)k * t * g1(x) - (k - 1) / (double)k * g2(x);
    };
    return func;
}

std::function<double(double)> LeastSquaresApproximator::orthogonalBasisFit(std::vector<std::function<double(double)>> funcs,
 std::function<double(double)> w){
    int length = funcs.size();
    std::vector<double> ks;

    for (int i = 0; i < length; ++i){
        double up = this->innerProduct(funcs[i], nullptr, w);
        double down = this->innerProduct(funcs[i], funcs[i], w);
        ks.push_back(up / down);
    }

    auto result_func = [funcs, ks](double x){
        double result = 0;
        for (int i = 0; i < funcs.size(); ++i)
            result += funcs[i](x) * ks[i];
        return result;
    };
    return result_func;
}

std::function<double(double)> LeastSquaresApproximator::normalEquationsFit(std::vector<std::function<double(double)>> funcs, 
std::function<double(double)> w){
    int length = funcs.size();
    Matrix A(length, length, 0);
    Matrix b(length, 1, 0);
    for (int i = 0; i < length; ++i){
        for (int j = i; j < length; ++j){
            A(i, j) = this->innerProduct(funcs[i], funcs[j], w);
            A(j, i) = A(i, j);
        }
        b(i, 0) = this->innerProduct(funcs[i], nullptr, w);
    }

    // 利用共轭梯度法求解正则方程
    IterationSolver solver(A, b);
    Matrix c = solver.ConjugateGradientSolve();

    auto result_func = [funcs, c](double x){
    double result = 0;
    for (int i = 0; i < funcs.size(); ++i)
        result += funcs[i](x) * c(i, 0);
    return result;
    };
    return result_func;
}

std::function<double(double)> LeastSquaresApproximator::LeastSquaresApproximator::standardPolynomial(unsigned int n){
    // 构造一组维度为n的k次多项式组成的基底
    std::vector<std::function<double(double)>> funcs;
    for (int i = 0; i < n + 1; ++i){
        auto func = this->standardPolynomialBase(i);
        funcs.push_back(func);
    }

    return this->normalEquationsFit(funcs);
}

std::function<double(double)> LeastSquaresApproximator::legendrePolynomial(unsigned int n){
    std::vector<std::function<double(double)>> funcs;
    // 构造勒让德多项式前两项
    auto f1 = [](double x){
        return 1;
    };
    auto f2 = [this](double x){
        double t = 2 * (x - this->min_x) / (this->max_x - this->min_x) - 1;
        return t;
    };
    funcs.push_back(f1), funcs.push_back(f2);

    for (int i = 2; i < n + 1; ++i){
        auto func = this->legendrePolynomialBase(i, funcs[i - 1], funcs[i - 2]);
        funcs.push_back(func);
    }

    return this->normalEquationsFit(funcs);
}

void Chapter05Demos::Demo1(){
    Matrix x = {{.1}, {.2}, {.3}, {.4}, {.5}, {.6}, {.7}, {.8}, {.9}};
    Matrix y = {{5.1234}, {5.3057}, {5.5687}, {5.9375}, {6.4370}, 
    {7.0978}, {7.9493}, {9.0253}, {10.3627}};

    LeastSquaresApproximator tool(x, y);

    // k次多项式拟合
    Matrix y1(x.rows, 1, 0);
    auto func1 = tool.standardPolynomial(4);
    std::cout << "k次多项式拟合结果" << std::endl;
    for (int i = 0; i < x.rows; ++i){
        y1(i, 0) = func1(x(i, 0));
        std::cout << y1(i, 0) << ", ";
    }
    std::cout << std::endl;
    std::cout << "误差为：" << norm(y - y1);
    std::cout << std::endl << std::endl;

    // 勒让德多项式拟合
    Matrix y2(x.rows, 1, 0);
    auto func2 = tool.legendrePolynomial(4);
    std::cout << "勒让德多项式拟合结果" << std::endl;
    for (int i = 0; i < x.rows; ++i){
        y2(i, 0) = func2(x(i, 0));
        std::cout << y2(i, 0) << ", ";
    }
    std::cout << std::endl;
    std::cout << "误差为：" << norm(y - y2);
    std::cout << std::endl << std::endl;
}
