# pragma once
#include <iostream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <iomanip>
#include <cmath>
#include <algorithm>

class Matrix{    
    public:
    size_t rows;  // 行数
    size_t cols;  // 列数

    // 矩阵的默认构造函数，构造未初始化的矩阵
    Matrix();

    // 矩阵构造函数
    Matrix(size_t rows, size_t cols, double default_value = 0.0);
    Matrix(std::initializer_list<std::initializer_list<double>> init);

    // 索引矩阵元素
    double& operator()(unsigned int r, unsigned int c);

    // 索引常量矩阵元素（不可修改矩阵）
    const double& operator() (unsigned int r, unsigned int c) const;

    // 矩阵和矩阵的加法运算
    Matrix operator+ (const Matrix& other) const;

    // 定义矩阵与常数的加法运算。注意模板函数的声明和实现必须在一起
    template <typename T>
    typename std::enable_if<
        // T 是算术类型（int/double/float 等标量）
        std::is_arithmetic<typename std::remove_cv_t<typename std::remove_reference<T>::type>>::value &&
        // T 移除引用和 const 后，不是 Matrix 类型（排除矩阵）
        !std::is_same<typename std::remove_cv_t<typename std::remove_reference<T>::type>, Matrix>::value,
        Matrix
    >::type
    operator+ (T&& x) const{
        Matrix result(rows, cols);
        for (int i = 0; i < data.size(); ++i)
            result.data[i] = this->data[i] + x;
        return result;
    }

    // 矩阵之间的减法运算
    Matrix operator- (const Matrix& other) const;

    // 矩阵与常数的减法运算
    template <typename T>
    typename std::enable_if<
        std::is_arithmetic<typename std::remove_cv_t<typename std::remove_reference<T>::type>>::value &&
        !std::is_same<typename std::remove_cv_t<typename std::remove_reference<T>::type>, Matrix>::value,
        Matrix
    >::type
    operator- (T&& x) const{
        Matrix result(rows, cols);
        for (int i = 0; i < data.size(); ++i)
            result.data[i] = this->data[i] - x;
        return result;
    }

    // 矩阵的一元-运算符，取反
    Matrix operator-() const;

    // 矩阵的流式输出
    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);

    // 矩阵间的乘法
    Matrix operator* (const Matrix& other) const;

    // 矩阵与常数的乘法
    template <typename T>
    typename std::enable_if<
        std::is_arithmetic<typename std::remove_cv_t<typename std::remove_reference<T>::type>>::value &&
        !std::is_same<typename std::remove_cv_t<typename std::remove_reference<T>::type>, Matrix>::value,
        Matrix
    >::type
    operator* (T&& x) const{
        Matrix result(rows, cols);
        for (int i = 0; i < data.size(); ++i)
            result.data[i] = this->data[i] * x;
        return result;
    }

    // 矩阵间的除法
    Matrix operator/ (const Matrix& other) const;

    // 矩阵与常数的除法
    template <typename T>
    typename std::enable_if<
        std::is_arithmetic<typename std::remove_cv_t<typename std::remove_reference<T>::type>>::value &&
        !std::is_same<typename std::remove_cv_t<typename std::remove_reference<T>::type>, Matrix>::value,
        Matrix
    >::type
    operator/ (T&& x) const{
        if (std::abs(x) <= 1e-12)  // 除数很小的时候警告
            std::cerr << "警告！除数过小！(<1e-12)" << std::endl;
        Matrix result(rows, cols);
        for (int i = 0; i < data.size(); ++i)
            result.data[i] = this->data[i] / x;
        return result;
    }

    friend Matrix operator/ (const double& x, const Matrix& mat);

    // 赋值运算符
    Matrix& operator= (const Matrix& other);
    // 矩阵初始化的赋值运算符
    Matrix& operator=(std::initializer_list<std::initializer_list<double>> init);

    // 矩阵的上三角部分
    Matrix upper_tri() const;

    // 矩阵的下三角部分
    Matrix lower_tri() const;

    // 矩阵的对角部分
    Matrix diag() const;

    // 矩阵的转置
    Matrix transpose() const;

    // 矩阵根据索引交换行或列
    Matrix sort(int* indexs, int axis) const;

    // 根据高斯约当消去法求矩阵的逆
    Matrix inv() const;

    private:
    std::vector<double> data;
};

// 矩阵与常数的右加、减、乘、除运算
Matrix operator- (const double& x, const Matrix& mat);
Matrix operator+ (const double& x, const Matrix& mat);
Matrix operator* (const double& x, const Matrix& mat);
Matrix operator/ (const double& x, const Matrix& mat);

/* 矩阵对象的工具函数 */
double norm(const Matrix& mat);  // 计算向量的模长
Matrix concat(const Matrix& mat_a, const Matrix& mat_b, int axis = 0);  // 沿行或者列拼接矩阵
Matrix generateRandomMatrix(size_t n, size_t m, double min_val = -10.0, double max_val = 10.0);  // 生成指定维度的非奇异随机方阵
Matrix generatePositiveDefiniteMatrix(size_t n);  // 生成随机正定矩阵，利用B^T * B为正定矩阵的性质(B可逆)
