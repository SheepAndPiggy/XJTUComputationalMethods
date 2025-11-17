#include "chapter04.hpp"
#include "chapter02.hpp"
#include "utils.hpp"
#include <functional>

InterpolateTool::InterpolateTool(Matrix x, Matrix y):
x(x), y(y), n(x.rows){
    if (x.rows != y.rows)
        throw std::invalid_argument("插值时插值节点坐标和值维度应相同！");
}

double InterpolateTool::_langrangeFunction(Matrix xs, Matrix ys, double x){
    double n = xs.rows;
    double result = 0;

    for (int i = 0; i < n; ++i){
        double xi = xs(i, 0);
        double lx = ys(i, 0);

        for (int j = 0; j < n; ++j){
            double xj = xs(j, 0);
            if (j != i)
                lx = lx * (x - xj) / (xi - xj);
        }
        result = result + lx;
    }
    return result;
}

double InterpolateTool::_error(Matrix xs, Matrix ys, double x, double (*function)(Matrix, Matrix, double)){
    Matrix x1(xs.rows - 1, 1, 0);
    Matrix y1(xs.rows - 1, 1, 0);
    Matrix x2(xs.rows - 1, 1, 0);
    Matrix y2(xs.rows - 1, 1, 0);

    for (int i = 0; i < xs.rows - 1; ++i){
        x1(i, 0) = xs(i, 0);
        x2(i, 0) = xs(i + 1, 0);
        y1(i, 0) = ys(i, 0);
        y2(i, 0) = ys(i + 1, 0);
    }

    double p1 = (*function)(x1, y1, x);
    double p2 = (*function)(x2, y2, x);

    return (p2 - p1) / (xs(xs.rows - 1, 0) - xs(0, 0)) * (x - xs(0, 0));
}

std::vector<SepInterpolateResult> InterpolateTool::_sep_function(int seg_num, double (*function)(Matrix, Matrix, double)){
    if (seg_num > (this->n - 1))
        throw std::invalid_argument("插值区间个数不可以大于n-1！");
    int seg_size = this->n / seg_num;
    std::vector<SepInterpolateResult> result_funcs;

    for (int i = 0; i < this->n; i+=seg_size){
        int end = i + seg_size;
        if (i + seg_size >= this->n)
            end = this->n - 1;

        Matrix xs(end - i + 1, 1, 0);
        Matrix ys(end - i + 1, 1, 0);
        for (int j = i; j <= end; ++j){
            xs(j - i, 0) = this->x(j, 0);
            ys(j - i, 0) = this->y(j, 0);
        }
        auto result_func = [function, xs, ys](double x){
            return (*function)(xs, ys, x);
        };
        auto error_func = [this, function, xs, ys](double x){
            return this->_error(xs, ys, x, function);
        };
        double range[2] = {this->x(i, 0), this->x(end, 0)};
        SepInterpolateResult result(range, error_func, result_func);
        result_funcs.push_back(result);
    }
    return result_funcs;
}

InterpolateResult InterpolateTool::_function(int seg_num, double (*function)(Matrix, Matrix, double)){
    std::vector<SepInterpolateResult> result_func = 
        this->_sep_function(seg_num, function);
    
    auto func = [result_func](double x){
        for (SepInterpolateResult f : result_func){
            if (f.in_range(x)){
                std::cout << "截断误差为：" << f.error(x) << std::endl;
                return f.interpolate(x);
            }
        }
        throw std::invalid_argument("数据不在插值区间范围内！");
    };

    auto error = [result_func](double x){
        for (SepInterpolateResult f : result_func){
            if (f.in_range(x))
                return f.error(x);
        }
        throw std::invalid_argument("数据不在插值区间范围内！");
    };

    InterpolateResult result(func, error);
    return result;
}

bool SepInterpolateResult::in_range(double x){
    if (x >= this->range[0] && x <= this->range[1])
        return true;
    return false;
}

InterpolateResult InterpolateTool::lagrangeInterpolation(int seg_num){
    InterpolateResult result = this->_function(seg_num, this->_langrangeFunction);
    return result;
}

double InterpolateTool::_newtonFunction(Matrix xs, Matrix ys, double x){
    int n = xs.rows;
    double* divided = new double[n];
    double* temp_divided = new double[n - 1]; // 初始化当前差商列
    double* last_temp_double = new double[n]; // 初始化上一个差商列
    for (int i = 0; i < n; ++i)
        last_temp_double[i] = ys(i, 0);  // 将上一个差商列赋值为f(x)
    divided[0] = last_temp_double[0];  // 将第一个牛顿多项式差商赋值为f(x0)

    for (int i = 1; i < n; ++i){        
        for (int j = 0; j < n - i; ++j)
            // 计算当前差商列
            temp_divided[j] = (last_temp_double[j + 1] - last_temp_double[j]) / (xs(j + i, 0)- xs(j, 0));
        divided[i] = temp_divided[0];  // 为牛顿多项式所需的差商赋值
        
        delete[] last_temp_double;
        last_temp_double = temp_divided;  // 将当前差商列转换为上一个差商列
        temp_divided = new double[n - i - 1];  // 重新初始化当前差商列
    }
    // 析构中间变量
    delete[] last_temp_double;
    delete[] temp_divided;

    // 计算插值结果
    double temp = 1;
    double result = 0;
    for (int i = 0; i < n; ++i){
        double f = divided[i];
        result += f * temp;
        temp *= x - xs(i, 0);
    }
    return result;
}

InterpolateResult InterpolateTool::newtonInterpolation(int seg_num){
    InterpolateResult result = this->_function(seg_num, this->_newtonFunction);
    return result;
}

InterpolateResult InterpolateTool::cubicSplineInterpolation(double* values){
    Matrix A(this->n - 2, this->n - 2, 0);
    Matrix b(this->n - 2, 1, 0);

    for (int i = 0; i < this->n - 2; ++i){
        double h1 = this->x(i + 1, 0) - this->x(i, 0);
        double h2 = this->x(i + 2, 0) - this->x(i + 1, 0);

        double mu = h1 / (h1 + h2);
        double lambda = 1 - mu;
        double d = 6 / (h1 + h2) * ((this->y(i + 2, 0) - this->y(i + 1, 0)) / h2 -
        (this->y(i + 1, 0) - this->y(i, 0)) / h1 );
        
        A(i, i) = 2;
        A(i, i == 0 ? i + 1 : i - 1) = i == 0 ? lambda : mu;
        A(i, i == this->n - 3 ? i - 1 : i + 1) = i == this->n - 3 ? mu : lambda;
        if (i == 0)
            b(i, 0) = d - mu * values[0];
        else if (i == this->n - 3)
            b(i, 0) = d - lambda * values[1];
        else
            b(i, 0) = d;
    }

    Matrix M = Matrix(this->n, 1, 0);
    M(0, 0) = values[0];
    M(this->n - 1, 0) = values[1];
    Matrix M_ = ThomasSolver::thomasSolve(A, b);
    for (int i = 1; i < this->n - 1; ++i)
        M(i, 0) = M_(i - 1, 0);
    
    auto func = [this, M](double x){
        for (int i = 0; i < this->n - 1; ++i){
            double x1 = this->x(i, 0);
            double x2 = this->x(i + 1, 0);
            double h = x2 - x1;
            double y1 = this->y(i, 0);
            double y2 = this->y(i + 1, 0);
            if (x >= x1 && x < x2){
                double m1 = M(i, 0);
                double m2 = M(i + 1, 0);
                double S = std::pow(x2 - x, 3) / 6 / h * m1 +
                std::pow(x - x1, 3) / 6 / h * m2 +
                (y1 - std::pow(h, 2) / 6 * m1) * (x2 - x) / h +
                (y2 - std::pow(h, 2) / 6 * m2) * (x - x1) / h;
                return S;
            }
        }
        throw std::invalid_argument("x不在插值允许区间内！");
    };

    auto error = [](double x){
        std::cout << "暂不计算三次样条插值的截断误差估计！" << std::endl;
        return 0;
    };

    InterpolateResult result(func, error);
    return result;
}

void Chapter04Demos::Demo1(int n){
    Matrix xs(n + 1, 1, 0);
    Matrix ys(n + 1, 1, 0);

    for (int i = 0; i < n + 1; ++i){
        xs(i, 0) = -1 + i / (double)n * 2;
        ys(i, 0) = 1 / (1 + 25 * std::pow(xs(i, 0), 2));
    }

    InterpolateTool tool(xs, ys);
    double test_xs[5] = {0.01, 0.13, 0.43, 0.73, 0.99};

    InterpolateResult lagrange_result = tool.lagrangeInterpolation();
    std::cout << std::string(60, '*') << std::endl;
    for (double x : test_xs){
        std::cout << "拉格朗日插值（不分段）在x=" << x << "处的插值结果："  << std::endl; 
        double ix = lagrange_result.interpolate(x);
        double iy = 1 / (1 + 25 * std::pow(x, 2));
        std::cout << "插值结果为：" << ix << std::endl;
        std::cout << "真实结果为：" << iy << std::endl;
        std::cout << "相对误差为：" << std::abs(ix - iy) / std::abs(iy) * 100 << "%" << std::endl << std::endl;
    }
    std::cout << std::string(60, '*') << std::endl << std::endl;

    InterpolateResult newton_result = tool.newtonInterpolation();
    std::cout << std::string(60, '*') << std::endl;
    for (double x : test_xs){
        std::cout << "牛顿插值（不分段）在x=" << x << "处的插值结果："  << std::endl; 
        double ix = newton_result.interpolate(x);
        double iy = 1 / (1 + 25 * std::pow(x, 2));
        std::cout << "插值结果为：" << ix << std::endl;
        std::cout << "真实结果为：" << iy << std::endl;
        std::cout << "相对误差为：" << std::abs(ix - iy) / std::abs(iy) * 100 << "%" << std::endl << std::endl;
    }
    std::cout << std::string(60, '*') << std::endl << std::endl;

    int sep_num = (int)(n / 3);
    InterpolateResult lagrange_result_sep = tool.lagrangeInterpolation(sep_num);
    std::cout << std::string(60, '*') << std::endl;
    for (double x : test_xs){
        std::cout << "拉格朗日插值（分" << sep_num << "段）在x=" << x << "处的插值结果："  << std::endl; 
        double ix = lagrange_result_sep.interpolate(x);
        double iy = 1 / (1 + 25 * std::pow(x, 2));
        std::cout << "插值结果为：" << ix << std::endl;
        std::cout << "真实结果为：" << iy << std::endl;
        std::cout << "相对误差为：" << std::abs(ix - iy) / std::abs(iy) * 100 << "%" << std::endl << std::endl;
    }
    std::cout << std::string(60, '*') << std::endl << std::endl;

    InterpolateResult newton_result_sep = tool.newtonInterpolation(sep_num);
    std::cout << std::string(60, '*') << std::endl;
    for (double x : test_xs){
        std::cout << "牛顿插值（分" << sep_num << "段）在x=" << x << "处的插值结果："  << std::endl; 
        double ix = lagrange_result_sep.interpolate(x);
        double iy = 1 / (1 + 25 * std::pow(x, 2));
        std::cout << "插值结果为：" << ix << std::endl;
        std::cout << "真实结果为：" << iy << std::endl;
        std::cout << "相对误差为：" << std::abs(ix - iy) / std::abs(iy) * 100 << "%" << std::endl << std::endl;
    }
    std::cout << std::string(60, '*') << std::endl << std::endl;

    double values[2] = {0, 0};
    InterpolateResult cubic_sep = tool.cubicSplineInterpolation(values);
    std::cout << std::string(60, '*') << std::endl;
    for (double x : test_xs){
        std::cout << "三次样条插值在x=" << x << "处的插值结果："  << std::endl; 
        double ix = cubic_sep.interpolate(x);
        double iy = 1 / (1 + 25 * std::pow(x, 2));
        std::cout << "插值结果为：" << ix << std::endl;
        std::cout << "真实结果为：" << iy << std::endl;
        std::cout << "相对误差为：" << std::abs(ix - iy) / std::abs(iy) * 100 << "%" << std::endl << std::endl;
    }
    std::cout << std::string(60, '*') << std::endl << std::endl;
}
