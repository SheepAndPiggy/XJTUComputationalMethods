#pragma once
#include <iostream>
#include "utils.hpp"
#include <cstring>
#include <chrono>
#include <cmath>

struct luResult{
    Matrix LU_mat;
    Matrix LU_P;

    luResult(Matrix& LU_mat, Matrix& LU_P);
    luResult(size_t rows, size_t cols);

    luResult& operator= (const luResult& other);
};

class DoolittleSolver{
    public: 
    // lu分解（列主元），返回包含L和U的矩阵和列主元的行索引
    static luResult luDecompose(const Matrix& mat, bool pivot = false);

    // 根据上或下三角矩阵进行回代，求解方程组的解，返回方程组的解
    static Matrix solveByTri(const Matrix& tri, const Matrix& b, std::string tri_type);
};

class SqrtMethodSolver{
    public:
    // 楚列斯基分解，返回分解后的矩阵G
    static Matrix choleskyDecompose(const Matrix& mat);

    // 改进平方根法，返回分解后的矩阵LU
    static Matrix improvedSqrtDecompose(const Matrix& mat);
};

class ThomasSolver{
    public:
    // 三对角追赶法
    static Matrix thomasSolve(const Matrix& mat, const Matrix& b);
};

struct QRResult{
    Matrix Q;
    Matrix R;

    QRResult(const Matrix& Q, const Matrix& R);
};

class RegularTransformer{
    public:
    // 构建维度为n，位置(i, j)的吉文斯变换矩阵
    static Matrix givensMatrix(size_t n, size_t i, size_t j, double s, double c);

    // 吉文斯变换的QR分解
    static QRResult givensQRDecompose(const Matrix& mat);

    // 豪斯霍尔德变换的QR分解
    static QRResult householderQRDecompose(const Matrix& mat);
};

// 测试用例
void Demo1();  // 计算实习2.1
void Demo11(size_t n = 1000);  // 测试LU分解分解正定矩阵和求解的性能
void Demo2();  // 计算实习2.2
void Demo21(size_t n = 1000);  // 比较LU分解、楚列斯基分解、改进平方根法的计算速度
void Demo3();  // 计算实习2.3
void Demo4();  // 计算实习2.4
void Demo41();  // 使用吉文斯变换进行QR分解，计算计算实习2.4