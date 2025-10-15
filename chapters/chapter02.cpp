#include "chapter02.hpp"
#include <iostream>
#include "utils.hpp"
#include <cstring>
#include <chrono>
#include <cmath>

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
            result(i, j) = mat(i, j) - temp;
            result(j, i) = result(i, j) / result(i, i);
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

QRResult::QRResult(const Matrix& Q, const Matrix& R): Q(Q), R(R){}

Matrix RegularTransformer::givensMatrix(size_t n, size_t i, size_t j, double s, double c){
    if (i >= j)
        throw std::invalid_argument("吉文斯变换参数i必须小于j！");
    Matrix givensMat(n, n, 0);
    for (int i_ = 0; i_ < n; ++i_)
        givensMat(i_, i_) = 1;
    givensMat(i, i) = c;
    givensMat(j, j) = c;
    givensMat(i, j) = s;
    givensMat(j, i) = -s;
    return givensMat;
}

 QRResult RegularTransformer::givensQRDecompose(const Matrix& mat){
    Matrix R = mat;
    Matrix Q(R.rows, R.rows, 0);  // 初始化Q矩阵
    for (int i = 0; i < R.rows; ++i)
        Q(i, i) = 1;

    for (int j = 0; j < R.cols; ++j){
        double first_sum = std::pow(R(j, j), 2);

        Matrix P_tot(R.rows - j, R.rows - j, 0);  // 初始化P矩阵
        for (int i = 0; i < R.rows - j; ++i)
            P_tot(i, i) = 1;

        for (int i = j + 1; i < R.rows; ++i){
            double first_sum_bak = first_sum;
            first_sum += std::pow(R(i, j), 2);

            double c = 0;
            if (i == j + 1)  // 第一次计算时c和xjj符号有关
                c = R(j, j) / std::sqrt(first_sum);
            else
                c = std::sqrt(first_sum_bak / first_sum);

            double s = R(i, j) / std::sqrt(first_sum);
            Matrix P = RegularTransformer::givensMatrix(R.rows - j, 0, i - j, s, c);
            P_tot = P * P_tot;
        }
        
        if (j != 0){
            Matrix I(j, j);
            for (int i = 0; i < j; ++i)
                I(i, i) = 1;
            Matrix Z(j, R.rows - j, 0);
            Matrix P_tot_concat = concat(concat(I, Z, 1), concat(Z.transpose(), P_tot, 1), 0);
            Q = P_tot_concat * Q;
            R = P_tot_concat * R;
        } else {
            Q = P_tot * Q;
            R = P_tot * R;
        }
    }

    QRResult result(Q.transpose(), R);
    return result;
 }

 QRResult RegularTransformer::householderQRDecompose(const Matrix& mat){
    Matrix R = mat;
    Matrix Q(R.rows, R.rows, 0);  // 初始化Q矩阵
    for (int i = 0; i < R.rows; ++i)
        Q(i, i) = 1;

    for (int j = 0; j < R.cols; ++j){
        double k = R(j, j) > 0 ? 1 : -1;
        Matrix line(R.rows - j, 1);  // 列向量
        for (int k = 0; k < R.rows - j; ++k)
            line(k, 0) = R(k + j, j); 
        double sigma = -k * norm(line);  // 计算sigma
        Matrix e(R.rows - j, 1);  // 单位向量
        e(0, 0) = 1;

        Matrix u = (line - sigma * e) / norm(line - sigma * e);  // 计算u
        Matrix I(R.rows - j, R.rows - j);  // 单位矩阵
        for (int k = 0; k < R.rows - j; ++k)
            I(k, k) = 1;
        Matrix H_ = I - 2 * u * u.transpose();

        if (j != 0){
            Matrix E(j, j);
            for (int k = 0; k < j; ++k)
                E(k, k) = 1;
            Matrix Z(j, R.rows - j, 0);
            Matrix H = concat(concat(E, Z, 1), concat(Z.transpose(), H_, 1), 0);

            Q = H * Q;
            R = H * R;
        } else{
            Q = H_ * Q;
            R = H_ * R;
        }
    }
    QRResult result(Q.transpose(), R);
    return result;
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
    
    luResult result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵

    std::cout << "LU分解矩阵：" << std::endl;
    std::cout << result.LU_mat << std::endl << std::endl;
    std::cout << "置换矩阵P：" << std::endl;
    std::cout << result.LU_P << std::endl;

    Matrix L = result.LU_mat.lower_tri();
    Matrix U = result.LU_mat.upper_tri();
    for (int i = 0; i < L.rows; ++i)
        L(i, i) = 1;  // 将L矩阵的对角线元素置为1

    Matrix y = DoolittleSolver::solveByTri(L, result.LU_P * b, "low");  // 求解Ly = Pb
    Matrix x = DoolittleSolver::solveByTri(U, y, "up");  // 求解Ux = y

    x = result.LU_P.transpose() * x;  // 根据行置换矩阵P交换结果x，置换矩阵的逆是其转置矩阵

    std::cout << "求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl << std::endl;
}

void Demo11(size_t n){
    std::cout << "矩阵维度：" << n << std::endl;

    Matrix A = generatePositiveDefiniteMatrix(n);

    Matrix b = generateRandomMatrix(n, 1);
    auto lu_start = std::chrono::high_resolution_clock::now();  // LU分解开始时间
    luResult result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵
    auto lu_end = std::chrono::high_resolution_clock::now();  // LU分解结束时间
    auto lu_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(lu_end - lu_start).count();

    std::cout << "LU分解用时（s）：" << std::endl;
    std::cout << lu_duration_ns / 1e9 << std::endl;  // 转换成秒

    Matrix L = result.LU_mat.lower_tri();
    Matrix U = result.LU_mat.upper_tri();
    for (int i = 0; i < L.rows; ++i)
        L(i, i) = 1;  // 将L矩阵的对角线元素置为1

    auto gauss_start = std::chrono::high_resolution_clock::now();  // 求解方程组开始时间
    Matrix y = DoolittleSolver::solveByTri(L, result.LU_P * b, "low");  // 求解Ly = Pb
    Matrix x = DoolittleSolver::solveByTri(U, y, "up");  // 求解Ux = y
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

    Matrix G = SqrtMethodSolver::choleskyDecompose(A);

    std::cout << "楚列斯基分解矩阵G：" << std::endl;
    std::cout << G << std::endl;

    Matrix LU = SqrtMethodSolver::improvedSqrtDecompose(A);

    std::cout << "改进平方根法分解矩阵LU：" << std::endl;
    std::cout << LU << std::endl;
}

void Demo21(size_t n){
        std::cout << "矩阵维度：" << n << std::endl;

    Matrix A = generatePositiveDefiniteMatrix(n);
    
    /* LU分解 */
    auto lu_start = std::chrono::high_resolution_clock::now();  // LU分解开始时间
    luResult lu_result = DoolittleSolver::luDecompose(A, true);  // 通过列主元的方法获取LU分解矩阵
    auto lu_end = std::chrono::high_resolution_clock::now();  // LU分解结束时间
    auto lu_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(lu_end - lu_start).count();

    std::cout << "LU分解用时（s）：" << std::endl;
    std::cout << lu_duration_ns / 1e9 << std::endl;  // 转换成秒

    /* 楚列斯基分解 */
    auto cholesky_start = std::chrono::high_resolution_clock::now();
    Matrix cholesky_result = SqrtMethodSolver::choleskyDecompose(A);
    auto cholesky_end = std::chrono::high_resolution_clock::now();
    auto cholesky_duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(cholesky_end - cholesky_start).count();

    std::cout << "楚列斯基分解用时（s）：" << std::endl;
    std::cout << cholesky_duration_ns / 1e9 << std::endl;

    /* 改进平方根法分解 */
    auto sqrt_start = std::chrono::high_resolution_clock::now();
    Matrix sqrt_result = SqrtMethodSolver::improvedSqrtDecompose(A);
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

    Matrix x = ThomasSolver::thomasSolve(A, b);

    std::cout << "三对角矩阵追赶法求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl;
}

void Demo4(){
    Matrix A(7, 7);
    Matrix b(7, 1);

    A = {
        {5, 4, 7, 5, 6, 7, 5}, 
        {4, 12, 8, 7, 8, 8, 6}, 
        {7, 8, 10, 9, 8, 7, 7}, 
        {5, 7, 9, 11, 9, 7, 5}, 
        {6, 8, 8, 9, 10, 8, 9}, 
        {7, 8, 7, 7, 8, 10, 10}, 
        {5, 6, 7, 5, 9, 10, 10}
    };
    b = {{39}, {53}, {56}, {53}, {58}, {57}, {52}};

    QRResult result = RegularTransformer::householderQRDecompose(A);
    Matrix y = result.Q.transpose() * b;
    Matrix x = DoolittleSolver::solveByTri(result.R, y, "up");

    std::cout << "利用豪斯霍尔德变换进行QR分解结果：" << std::endl;
    std::cout << "Q矩阵：" << std::endl;
    std::cout << result.Q << std::endl;
    std::cout << "R矩阵：" << std::endl;
    std::cout << result.R << std::endl;

    std::cout << "求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl;
}

void Demo41(){
    Matrix A(7, 7);
    Matrix b(7, 1);

    A = {
        {5, 4, 7, 5, 6, 7, 5}, 
        {4, 12, 8, 7, 8, 8, 6}, 
        {7, 8, 10, 9, 8, 7, 7}, 
        {5, 7, 9, 11, 9, 7, 5}, 
        {6, 8, 8, 9, 10, 8, 9}, 
        {7, 8, 7, 7, 8, 10, 10}, 
        {5, 6, 7, 5, 9, 10, 10}
    };
    b = {{39}, {53}, {56}, {53}, {58}, {57}, {52}};

    QRResult result = RegularTransformer::givensQRDecompose(A);
    Matrix y = result.Q.transpose() * b;
    Matrix x = DoolittleSolver::solveByTri(result.R, y, "up");

    std::cout << "利用吉文斯变换进行QR分解结果：" << std::endl;
    std::cout << "Q矩阵：" << std::endl;
    std::cout << result.Q << std::endl;
    std::cout << "R矩阵：" << std::endl;
    std::cout << result.R << std::endl;

    std::cout << "求解结果：" << std::endl;
    std::cout << x.transpose() << std::endl;
}
