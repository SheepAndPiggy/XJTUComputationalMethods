#include "utils.hpp"
#include "chapter03.hpp"

IterationSolver::IterationSolver(const Matrix& A, const Matrix& b, unsigned int max_iter, double epsilon):
 A(A), b(b), max_iter(max_iter), epsilon(epsilon){
    this->D = A.diag();

    this->D_inv = this->D;
    for (int i = 0; i < A.rows; ++i)
        this->D_inv(i, i) = 1 / this->D(i, i);
    
    this->E = -(A.lower_tri() - A.diag());
    this->F = -(A.upper_tri() - A.diag());
}

Matrix IterationSolver::baseSolve(const Matrix& B, const Matrix& g){
    // 使用0-1的随机数进行初始化
    Matrix xk = generateRandomMatrix(B.cols, 1, 1e-10, 1);
    for (int i = 0; i < max_iter; ++i){
        Matrix xk_old = xk;
        xk = B * xk + g;
        double error = this->convergenceError(xk_old, xk);
        double final_error = this->residualNormError(xk) / norm(xk);
        std::cout << "第" << i + 1 << "步" << "残差||Ax-b||/||x||：" << std::endl;
        std::cout << final_error << std::endl;
        if (error <= epsilon){
            std::cout << "*****求解收敛！*****" << std::endl;
            if (final_error <= epsilon)
                std::cout << "*****求解成功！*****" << std::endl;
            else
                std::cout << "*****收敛但求解不成功！*****" << std::endl;
            return xk;
        }
    }
    std::cout << "*****求解不收敛！*****" << std::endl;
    return xk;
}

Matrix IterationSolver::jacobiSolve(){
    Matrix B = D_inv * (E + F);
    Matrix g = D_inv * b;
    Matrix x = this->baseSolve(B, g);
    return x;
}

Matrix IterationSolver::SORSolve(double w){
    Matrix B = (D - w * E).inv() * ((1 - w) * D + w * F);
    Matrix g = w * (D - w * E).inv() * b;
    Matrix x = this->baseSolve(B, g);
    return x;
}

Matrix IterationSolver::gaussSolve(){
    return this->SORSolve(1);
}

double IterationSolver::convergenceError(const Matrix& x_old, const Matrix& x_new){
    double error = norm((x_new - x_old) / x_old);  // 相对误差
    return error;
}

double IterationSolver::residualNormError(const Matrix& x){
    double error = norm(A * x - b);
    return error;
}

void Chapter03Demos::Demo1(){
    Matrix A = {
        {2.52, 0.95, 1.25, -0.85}, 
        {0.39, 1.69, -0.45, 0.49}, 
        {0.55, -1.25, 1.96, -0.98}, 
        {0.23, -1.15, -0.45, 2.31}
    };
    Matrix b = {{1.38}, {-0.34}, {0.67}, {1.52}};

    IterationSolver solver(A, b, 10000, 1e-3);

    std::cout << solver.jacobiSolve().transpose() << std::endl;;
    std::cout << solver.gaussSolve().transpose() << std::endl;
}