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

    // 利用二分法找到解存在区间
    std::vector<std::pair<double, double>> findAllRootIntervals(
    double L, double R,
    double step);

    private:
    std::function<double(double)> func;
};

struct Chapter07Demos{
    static void Demo1();
};
