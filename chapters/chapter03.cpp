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

Matrix IterationSolver::ConjugateGradientSolve(){
    Matrix x = generateRandomMatrix(A.cols, 1, -1, 1);

    Matrix r = b - A * x;
    Matrix d = r;
    double beta;

    for (int i = 0; i < max_iter; ++i){
        double final_error = this->residualNormError(x) / norm(x);
        // std::cout << "第" << i + 1 << "步" << "残差||Ax-b||/||x||：" << std::endl;
        // std::cout << final_error << std::endl;
        if (final_error < epsilon){
            std::cout << "维度为" << A.cols << "的共轭梯度法过程共用" << i + 1 << "步, "
             << "残差||Ax-b||/||x||为" << final_error << std::endl;
            break;
        }
        double alpha = (r.transpose() * d)(0, 0) / 
        (d.transpose() * A * d)(0, 0);
        x = x + alpha * d;
        r = b - A * x;
        beta = -(r.transpose() * A * d)(0, 0) / 
        (d.transpose() * A * d)(0, 0);
        d = r + beta * d;
    }
    return x;
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

void Chapter03Demos::Demo2(){
    int ns[3] = {100, 200, 400};

    for (int n: ns){
        Matrix A(n, n, 0);
        Matrix b(n, 1, 0);
        for (int i = 0; i < n; ++i){
            int j_start = (i > 0) ? (i - 1) : 0;
            int j_end   = (i < n - 1) ? (i + 1) : i;
            for (int j = j_start; j <= j_end; ++j)
                A(i, j) = (i == j) ? -2 : 1;
        }
        b(0, 0) = -1;
        b(n - 1, 0) = -1;
        
        IterationSolver solver(A, b, 10000, 1e-3);
        Matrix result = solver.ConjugateGradientSolve();
    }
}