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
        k +=1;
    }
    std::cout << "迭代法共用" << k << "步" << std::endl;
    return x2;
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
}
