#pragma once
#include "utils.hpp"
#include "chapter02.hpp"
#include <functional>

struct SepInterpolateResult{
    double range[2];
    std::function<double(double)> error;  // 截断区间误差计算函数
    std::function<double(double)> interpolate;  // 区间插值函数

    SepInterpolateResult(double range_[2], std::function<double(double)> error, std::function<double(double)> interpolate): 
    error(error), interpolate(interpolate){
        range[0] = range_[0];
        range[1] = range_[1];
    };

    bool in_range(double x);  // 判断x是否在区间插值函数的范围内
};

struct InterpolateResult{
    std::function<double(double)> error;  // 截断误差计算函数
    std::function<double(double)> interpolate;  // 插值函数

    InterpolateResult(std::function<double(double)> interpolate, std::function<double(double)> error):
    error(error), interpolate(interpolate){}
};

class InterpolateTool{
    public:
    InterpolateTool(Matrix x, Matrix y);

    InterpolateResult lagrangeInterpolation(int seg_num = 1);
    InterpolateResult newtonInterpolation(int seg_num = 1);
    InterpolateResult cubicSplineInterpolation(double* values);

    private:
    double n;
    Matrix x;
    Matrix y;

    static double _langrangeFunction(Matrix xs, Matrix ys, double x);
    static double _newtonFunction(Matrix xs, Matrix ys, double x);

    std::vector<SepInterpolateResult> _sep_function(int seg_num, double (*function)(Matrix, Matrix, double));
    InterpolateResult _function(int seg_num, double (*function)(Matrix, Matrix, double));
    double _error(Matrix xs, Matrix ys, double x, double (*function)(Matrix, Matrix, double));
};

// 测试用例
struct Chapter04Demos{
    static void Demo1(int n);
};