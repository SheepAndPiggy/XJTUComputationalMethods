#include <iostream>
#include <cstdlib>  // 引用C标准库
#include <vector>  // 引用向量标准库
#include <type_traits>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "utils.hpp"

/* 矩阵（一维vector模拟二维矩阵） */
Matrix::Matrix(): rows(-1), cols(-1), data(std::vector<double>()){}

Matrix::Matrix(size_t rows, size_t cols, double default_value):  // 默认参数只能在类定义的时候传入一次
    rows(rows),
    cols(cols),
    data(rows * cols, default_value)
    {}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> init):
rows(-1), cols(-1), data(std::vector<double>()){
    int i = 0;
    int j = 0;
    for (const auto& row : init) {  // 遍历外层初始化列表（行）
        j = 0;
        for (double val : row) {  // 遍历内层初始化列表（列）
            this->data.push_back(val);  // 赋值到 data
            j++;
        }
        i++;
    }
    this->rows = i;
    this->cols = j;
}

double& Matrix::operator()(unsigned int r, unsigned int c){
    if ((r + 1) > rows || (c + 1) > cols)
        throw std::invalid_argument("索引超过矩阵范围！");
    return data[r * cols + c];
}

// 常量矩阵获取数值，'{'之前的const表示该函数不会修改对象内容，相当于隐式传入了const Matrix* this，而不是Matrix* this
const double& Matrix::operator() (unsigned int r, unsigned int c) const{
    if ((r + 1) > rows || (c + 1) > cols)
        throw std::invalid_argument("索引超过矩阵范围！");
    return data[r * cols + c];
}

// 矩阵的+运算符，时间复杂度为O(n)
Matrix Matrix::operator+ (const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("应用矩阵加法时，矩阵维度必须一致！");
        
    Matrix result(rows, cols);
    for (int i = 0; i < data.size(); ++i)
        result.data[i] = this->data[i] + other.data[i];
    return result;
}

Matrix Matrix::operator- (const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("应用矩阵减法时，矩阵维度必须一致！");
        
    Matrix result(rows, cols);
    for (int i = 0; i < data.size(); ++i)
        result.data[i] = this->data[i] - other.data[i];
    return result;
}

Matrix Matrix::operator-() const{
    Matrix result(rows, cols);
    for (int i = 0; i < data.size(); ++i)
        result.data[i] = -this->data[i];
    return result;
}

Matrix operator- (const double& x, const Matrix& mat){
    return -(mat - x);
}

Matrix operator+ (const double& x, const Matrix& mat){
    return mat + x;
}

Matrix operator* (const double& x, const Matrix& mat){
    return mat * x;
}

Matrix operator/ (const double& x, const Matrix& mat){
    Matrix result(mat.rows, mat.cols);
    for (int i = 0; i < result.data.size(); ++i){
        if (std::abs(mat.data[i]) <= 1e-12)
            std::cerr << "警告！除数矩阵存在接近0的元素！（<1e-12）" << std::endl;
        result.data[i] = x / mat.data[i];
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& mat) {
    for (size_t i = 0; i < mat.rows; ++i) {
        for (size_t j = 0; j < mat.cols; ++j) {
            os << std::setw(8) << std::fixed << std::setprecision(3)  // 宽度为8，3位小数
               << mat.data[i * mat.cols + j];
        }
        os << std::endl;  // 每一行结束后换行
    }
    return os;
}

Matrix Matrix::operator* (const Matrix& other) const{
    if (cols != other.rows)
        throw std::invalid_argument("矩阵间乘法维度不匹配！");
    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i){
        for (int k = 0; k < other.cols; ++k){
            double sum = 0.0;
            for (int j = 0; j < cols; ++j)
                sum += this->data[i * cols + j] * other.data[j * other.cols + k];
            result.data[i * result.cols + k] = sum;
        }
    }
    return result;
}

Matrix Matrix::operator/ (const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("矩阵间除法维度不匹配！");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j){
            if (std::abs(other(i, j)) < 1e-12)
                throw std::domain_error("矩阵含有零元素，无法作为除数！");
            result(i, j) = (*this)(i, j) / other(i, j);
        }
    return result;
}

Matrix& Matrix::operator= (const Matrix& other){
    if (this == &other)
        return *this;
    if (this->rows == -1 && this->cols == -1){  // 如果矩阵没有经过初始化
        this->rows = other.rows;
        this->cols = other.cols;
        this->data = other.data;  // 直接复制所有的值
        return *this;
    }
    if (this->rows != other.rows || this->cols != other.cols)
        throw std::invalid_argument("赋值符号'='两侧矩阵大小不一致！");
    std::copy(other.data.begin(), other.data.end(), this->data.begin());
    return *this;
}

Matrix& Matrix::operator=(std::initializer_list<std::initializer_list<double>> init){
    if (init.size() != rows)
        throw std::invalid_argument("初始化列表行数与矩阵行数不匹配");
    int i = 0;
    for (const auto& row : init) {  // 遍历外层初始化列表（行）
        if (row.size() != cols) {
            throw std::invalid_argument("初始化列表列数与矩阵列数不匹配");
        }

        int j = 0;
        for (double val : row) {  // 遍历内层初始化列表（列）
            this->data[i * cols + j] = val;  // 赋值到 data
            j++;
        }
        i++;
    }
    return *this;
}

Matrix Matrix::upper_tri() const{
    if (rows != cols)
        throw std::invalid_argument("非方阵的矩阵不存在上三角矩阵！");
    Matrix result(rows, cols, 0);
    for (int i = 0; i < rows; ++i)
        for(int j = i; j < cols; ++j)
            result.data[i * cols + j] = this->data[i * cols + j];
    return result;
}

Matrix Matrix::lower_tri() const{
    if (rows != cols)
        throw std::invalid_argument("非方阵的矩阵不存在下三角矩阵！");
    Matrix result(rows, cols, 0);
    for (int i = 0; i < rows; ++i)
        for(int j = 0; j <= i; ++j)
            result.data[i * cols + j] = this->data[i * cols + j];
    return result;
}

Matrix Matrix::diag() const{
    if (rows != cols)
        throw std::invalid_argument("非方阵的矩阵不存在对角矩阵！");
    Matrix result(rows, cols, 0);
    for (int i = 0; i < rows; ++i)
        result.data[i * cols + i] = this->data[i * cols + i];
    return result;
}

Matrix Matrix::transpose() const{
    Matrix result(cols, rows, 0);
    for (int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            result.data[j * rows + i] = this->data[i * cols + j];
    return result;
}

Matrix Matrix::sort(int* indexs, int axis) const{
    Matrix result(rows, cols);
    result = *this;
    if (axis == 0){  // 交换行
        for (int i = 0; i < rows; ++i){
            int index = indexs[i];
            if (index != i)
                for (int j = 0; j < cols; ++j)
                    result(index, j) = (*this)(i, j);
        }
    } else if (axis = 1){  // 交换列
        for (int i = 0; i < cols; ++i){
            int index = indexs[i];
            if (index != i)
                for (int j = 0; j < rows; ++j)
                    result(j, index) = (*this)(j, i);
        }
    } else{
        throw std::invalid_argument("axis参数只能为0(行方向)和1(列方向)！");
    }
    return result;
}

Matrix Matrix::inv() const{
    if (rows != cols)
        throw std::invalid_argument("求逆的矩阵必须为方阵！");
    
    Matrix I(rows, cols);
    for (int i = 0; i < rows; ++i)
        I(i, i) = 1;
    
    Matrix GY_mat = concat(*this, I, 1);

    for (int i = 0; i < rows; ++i){
        double q = GY_mat(i, i);
        for (int k = i; k < GY_mat.cols; ++k)
            GY_mat(i, k) /= q;

        for (int j = 0; j < rows; ++j){
            if (j == i)
                continue;
            double l = -GY_mat(j, i);
            for (int k = i; k < GY_mat.cols; ++k)
                GY_mat(j, k) += l * GY_mat(i, k);
        }
    }

    Matrix inv_mat(rows, cols);
    for (int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            inv_mat(i, j) = GY_mat(i, j + cols);
    return inv_mat;
}

/* 矩阵对象的工具函数 */
double norm(const Matrix& mat){
    double sum = 0;
    if (mat.rows == 1){
        for (int i = 0; i < mat.cols; ++i)
            sum += std::pow(mat(0, i), 2);
    } else if (mat.cols == 1){
        for (int i = 0; i < mat.rows; ++i)
            sum += std::pow(mat(i, 0), 2);
    } else {
        throw std::invalid_argument("求模长输入必须为向量！");
    }
    return std::sqrt(sum);
}

Matrix concat(const Matrix& mat_a, const Matrix& mat_b, int axis){
    size_t rows, cols;
    if (axis == 0){
        if (mat_a.cols != mat_b.cols)
            throw std::invalid_argument("垂直拼接的两个矩阵列方向维度不一致！");
        rows = mat_a.rows + mat_b.rows;
        cols = mat_a.cols;
        Matrix result(rows, cols);
        for (int j = 0; j < cols; ++j){
            for (int i = 0; i < mat_a.rows; ++i)
                result(i, j) = mat_a(i, j);
            for (int i = mat_a.rows; i < rows; ++i)
                result(i, j) = mat_b(i - mat_a.rows, j);
        }
        return result;
    } else if (axis == 1){
        if (mat_a.rows != mat_b.rows)
            throw std::invalid_argument("水平拼接的两个矩阵行方向维度不一致！");
        rows = mat_a.rows;
        cols = mat_a.cols + mat_b.cols;
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i){
            for (int j = 0; j < mat_a.cols; ++j)
                result(i, j) = mat_a(i, j);
            for (int j = mat_a.cols; j < cols; ++j)
                result(i, j) = mat_b(i, j - mat_a.cols);
        }
        return result;
    } else {
        throw std::invalid_argument("axis参数只能为0(行方向)和1(列方向)！");
    }
}

Matrix generateRandomMatrix(size_t n, size_t m, double min_val, double max_val) {
    Matrix B(n, m);
    std::srand(std::time(nullptr));  // 随机数种子

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            // 生成 [min_val, max_val) 之间的随机数
            double rand_val = min_val + (max_val - min_val) * std::rand() / RAND_MAX;
            B(i, j) = rand_val;
        }
    }

    int k = m;
    if (n <= m)  // 如果矩阵的行大于列
        k = n;

    for (size_t i = 0; i < k; ++i) {
        B(i, i) += n * std::abs(max_val);  // 对角线元素增强，降低奇异概率
    }
    return B;
}

Matrix generatePositiveDefiniteMatrix(size_t n) {
    Matrix B = generateRandomMatrix(n, n);  // 生成随机可逆矩阵 B
    return B.transpose() * B;              // A = B^T * B 是正定矩阵
}
