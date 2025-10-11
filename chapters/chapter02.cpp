#include "chapter02.hpp"
#include <iostream>
#include "utils.hpp"
#include <cstring>
#include <chrono>

luResult::luResult(Matrix& LU_mat, Matrix& LU_P):
LU_mat(LU_mat), LU_P(LU_P){}

luResult::luResult(size_t rows, size_t cols):
LU_mat(rows, cols), LU_P(rows, cols){}

luResult& luResult::operator= (const luResult& other){
    if (this == &other)
        return (*this);
    
    this->LU_mat = other.LU_mat;
    this->LU_P = other.LU_P;
    return (*this);
}

luResult DoolittleSolver::luDecompose(const Matrix& mat, bool pivot){
    Matrix result = mat;
    // 定义行索引列表
    int indexs[mat.rows];
    // 初始化行索引列表
    for (int i_ = 0; i_ < mat.rows; ++i_)
        indexs[i_] = i_;

    for (int i = 0; i < mat.rows - 1; ++i){
        if (pivot){
            // 找到主元所在行
            unsigned int max_index = i;
            double temp = std::abs(result(i, i));
            for (int i_ = i; i_ < mat.rows; ++i_){
                if (std::abs(result(i_, i)) > temp){
                    max_index = i_;
                    temp = std::abs(result(i_, i));
                }
            }

            // 将主元所在行与当前行进行交换，并保存行索引列表
            int temp_int = indexs[i];
            indexs[i] = indexs[max_index];
            indexs[max_index] = temp_int;
            for (int j_ = 0; j_ < mat.cols; ++j_){
                double temp_double = result(i, j_);
                result(i, j_) = result(max_index, j_);
                result(max_index, j_) = temp_double;
            }
        }
        
        // 判断主元是否接近0
        double main_value = result(i, i);
        if (std::abs(main_value) < 1e-12)
            throw std::domain_error("主元接近0！矩阵可能奇异！");

        for (int j = i + 1; j < mat.rows; ++j){
            double l = result(j, i) / main_value;
            for (int k = i; k < mat.cols; ++k){
                result(j, k) -= l * result(i, k);  // 将主元所在列之后消去为0
            }
            result(j, i) = l;  // 保存l_ij至矩阵下三角
        }
    }
    Matrix LU_P(mat.rows, mat.rows);
    for (int i = 0; i < mat.rows; ++i)
        LU_P(i, indexs[i]) = 1;
    luResult result_stract(result, LU_P);
    return result_stract;
}

Matrix DoolittleSolver::solveByTri(const Matrix& tri, const Matrix& b, std::string tri_type){
    if (tri.rows != tri.cols)
        throw std::invalid_argument("上/下三角矩阵回代求解，矩阵必须为方阵！");
    Matrix result(tri.rows, b.cols);
    if (tri_type == "up"){
        for (int k = 0; k < b.cols; ++k){  // 遍历b的每列
            for (int i = tri.rows - 1; i >= 0; --i){
                double temp = 0;
                for (int j = i + 1; j < tri.cols; ++j)
                    temp += tri(i, j) * result(j, k);
                result(i, k) = (b(i, k) - temp) / tri(i, i);
            }
        }
    } else if (tri_type == "low"){
        for (int k = 0; k < b.cols; ++k){
            for (int i = 0; i < tri.rows; ++i){
                double temp = 0;
                for (int j = 0; j < i; ++j)
                    temp += tri(i, j) * result(j, k);
                result(i, k) = (b(i, k) - temp) / tri(i, i);
            }
        }
    } else{
        throw std::invalid_argument("三角矩阵类型必须为'up'或者'low'！");
    }
    return result;
}

Matrix SqrtMethodSolver::choleskyDecompose(const Matrix& mat){
    if (mat.rows != mat.cols)
        throw std::invalid_argument("楚列斯基分解必须是方阵！");
    
    Matrix G(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; ++i){
        for (int j = 0; j < i; ++j){
            double temp = 0;
            for (int k = 0; k < j; ++k)
                temp += G(i, k) * G(j, k);
            G(i, j) = (mat(i, j) - temp) / G(j, j);
        }
        double temp = 0;
        for (int k = 0; k < i; ++k)
            temp += G(i, k) * G(i, k);
        G(i, i) = std::sqrt(mat(i, i) - temp);
    }
    return G;
}

Matrix SqrtMethodSolver::improvedSqrtDecompose(const Matrix& mat){
    if (mat.rows != mat.cols)
        throw std::invalid_argument("改进平方根法分解LDL^T必须是方阵！");
    
    Matrix result(mat.rows, mat.cols);
    for (int i = 0; i < mat.rows; ++i){
        double temp = 0;
        for (int k = 0; k < i; ++k)
            temp += result(i, k) * result(k, i);
        result(i, i) = mat(i, i) - temp;
        for (int j = i + 1; j < mat.cols; ++j){
            double temp = 0;
            for (int k = 0; k < i; ++k){
                temp += result(i, k) * result(k, j);
            }
            result(j, i) = (mat(i, j) - temp) / result(i, i);
        }
    }
    return result;
}

Matrix ThomasSolver::thomasSolve(const Matrix& mat, const Matrix& b){
    Matrix y(mat.rows, 1);
    Matrix u(mat.rows, 1);
    Matrix x(mat.rows, 1);

    u(0, 0) = mat(0, 0);
    y(0, 0) = b(0, 0);
    for (int i = 1; i < mat.rows; ++i){
        double l = mat(i, i - 1) / u(i - 1, 0);
        u(i, 0) = mat(i, i) - l * mat(i - 1, i);
        y(i, 0) = b(i, 0) - l * y(i - 1, 0);
    }

    x(mat.rows - 1, 0) = y(mat.rows - 1, 0) / u(mat.rows - 1, 0);
    for (int i = mat.rows - 2; i >= 0; --i)
        x(i, 0) = (y(i, 0) - mat(i, i + 1) * x(i + 1, 0)) / u(i, 0);
    return x;
}


void Demo1(){
    Matrix A(4, 4);
    A = {
        {1.1348, 3.8326, 1.1651, 3.4017}, 
        {0.5301, 1.7875, 2.5330, 1.5435}, 
        {3.4129, 4.9317, 8.7643, 1.3142}, 
        {1.2371, 4.9998, 10.6721, 0.0147}
    };

    Matrix b(4, 1);
    b = {
        {9.5342}, {6.3941}, {18.4231}, {16.9237}
    };
    
    luResult result(A.rows, A.cols);
    result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵

    std::cout << "LU分解矩阵：" << std::endl;
    std::cout << result.LU_mat << std::endl << std::endl;
    std::cout << "置换矩阵P：" << std::endl;
    std::cout << result.LU_P << std::endl;

    Matrix L(result.LU_mat.rows, result.LU_mat.cols);
    Matrix U(result.LU_mat.rows, result.LU_mat.cols);
    L = result.LU_mat.lower_tri();
    U = result.LU_mat.upper_tri();
    for (int i = 0; i < L.rows; ++i)
        L(i, i) = 1;  // 将L矩阵的对角线元素置为1
    
    Matrix y(result.LU_mat.cols, b.cols);
    Matrix x(result.LU_mat.cols, b.cols);

    y = DoolittleSolver::solveByTri(L, result.LU_P * b, "low");  // 求解Ly = Pb
    x = DoolittleSolver::solveByTri(U, y, "up");  // 求解Ux = y

    x = result.LU_P.transpose() * x;  // 根据行置换矩阵P交换结果x，置换矩阵的逆是其转置矩阵

    std::cout << "求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl << std::endl;
}

void Demo11(size_t n){
    std::cout << "矩阵维度：" << n << std::endl;

    Matrix A(n, n);
    A = generatePositiveDefiniteMatrix(n);

    Matrix b(n, 1);
    b = generateRandomMatrix(n, 1);
    
    luResult result(A.rows, A.cols);
    auto lu_start = std::chrono::high_resolution_clock::now();  // LU分解开始时间
    result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵
    auto lu_end = std::chrono::high_resolution_clock::now();  // LU分解结束时间
    auto lu_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(lu_end - lu_start).count();

    std::cout << "LU分解用时（s）：" << std::endl;
    std::cout << lu_duration_ns / 1e9 << std::endl;  // 转换成秒

    Matrix L(result.LU_mat.rows, result.LU_mat.cols);
    Matrix U(result.LU_mat.rows, result.LU_mat.cols);
    L = result.LU_mat.lower_tri();
    U = result.LU_mat.upper_tri();
    for (int i = 0; i < L.rows; ++i)
        L(i, i) = 1;  // 将L矩阵的对角线元素置为1
    
    Matrix y(result.LU_mat.cols, b.cols);
    Matrix x(result.LU_mat.cols, b.cols);

    auto gauss_start = std::chrono::high_resolution_clock::now();  // 求解方程组开始时间
    y = DoolittleSolver::solveByTri(L, result.LU_P * b, "low");  // 求解Ly = Pb
    x = DoolittleSolver::solveByTri(U, y, "up");  // 求解Ux = y
    auto gauss_end = std::chrono::high_resolution_clock::now();  // 求解方程组结束时间
    auto gauss_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(gauss_end - gauss_start).count();
    std::cout << "LU分解法方程组求解用时（s）：" << std::endl;
    std::cout << gauss_duration_ns / 1e9 << std::endl;  // 转换成秒
    std::cout << std::endl;
}

void Demo2(){
    Matrix A(20, 20);
    for (int i = 0; i < A.rows; ++i){  // 初始化矩阵A
        for (int j = 0; j < A.cols; ++j)
            A(i, j) = (i == j) ? i + 1 : std::min(i, j) + 1;
    }

    Matrix G(A.rows, A.cols);
    G = SqrtMethodSolver::choleskyDecompose(A);

    std::cout << "楚列斯基分解矩阵G：" << std::endl;
    std::cout << G << std::endl;

    Matrix LD(A.rows, A.cols);
    LD = SqrtMethodSolver::improvedSqrtDecompose(A);

    std::cout << "改进平方根法分解矩阵LU：" << std::endl;
    std::cout << LD << std::endl;
}

void Demo21(size_t n){
        std::cout << "矩阵维度：" << n << std::endl;

    Matrix A(n, n);
    A = generatePositiveDefiniteMatrix(n);
    
    /* LU分解 */
    luResult lu_result(A.rows, A.cols);
    auto lu_start = std::chrono::high_resolution_clock::now();  // LU分解开始时间
    lu_result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵
    auto lu_end = std::chrono::high_resolution_clock::now();  // LU分解结束时间
    auto lu_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(lu_end - lu_start).count();

    std::cout << "LU分解用时（s）：" << std::endl;
    std::cout << lu_duration_ns / 1e9 << std::endl;  // 转换成秒

    /* 楚列斯基分解 */
    Matrix cholesky_result(A.rows, A.cols);
    auto cholesky_start = std::chrono::high_resolution_clock::now();
    cholesky_result = SqrtMethodSolver::choleskyDecompose(A);
    auto cholesky_end = std::chrono::high_resolution_clock::now();
    auto cholesky_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(cholesky_end - cholesky_start).count();

    std::cout << "楚列斯基分解用时（s）：" << std::endl;
    std::cout << cholesky_duration_ns / 1e9 << std::endl;

    /* 改进平方根法分解 */
    Matrix sqrt_result(A.rows, A.cols);
    auto sqrt_start = std::chrono::high_resolution_clock::now();
    sqrt_result = SqrtMethodSolver::improvedSqrtDecompose(A);
    auto sqrt_end = std::chrono::high_resolution_clock::now();
    auto sqrt_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(sqrt_end - sqrt_start).count();

    std::cout << "改进平方根法分解用时（s）：" << std::endl;
    std::cout << sqrt_duration_ns / 1e9 << std::endl;
}

void Demo3(){
    Matrix A(20, 20, 0);
    Matrix b(20, 1, 2);

    b(0, 0) = 3;
    b(b.rows - 1, 0) = 3;
    for (int i = 0; i < A.rows; ++i){
        A(i, i) = 4;
        if (i > 0)
            A(i, i - 1) = -1;
        if (i < A.rows - 1)
            A(i, i + 1) = -1;
    }

    Matrix x(A.cols, 1);
    x = ThomasSolver::thomasSolve(A, b);

    std::cout << "三对角矩阵追赶法求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl;
}
