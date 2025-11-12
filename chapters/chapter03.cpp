#include "utils.hpp"
#include "chapter03.hpp"
#include "chapter02.hpp"

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

ArnoldResult IterationSolver::ArnoldiProcess(const Matrix& r, int m){
    Matrix H(m + 1, m, 0);
    Matrix* v = new Matrix[m + 1];
    
    v[0] = r / norm(r);  // 初始化v0

    // 构造上海森伯格矩阵Hm
    for (int j = 0; j < m; ++j){
        for (int i = 0; i <= j; ++i)
            H(i, j) = (v[i].transpose() * A * v[j])(0, 0);
        
        if (j < m){
            v[j + 1] = A * v[j];
            for (int k = 0; k <= j; ++k)
                v[j + 1] = v[j + 1] - H(k, j) * v[k];

            H(j + 1, j) = norm(v[j + 1]);
            v[j + 1] = v[j + 1] / H(j + 1, j);
        }
    }
    ArnoldResult result(v, H);
    return result;
}

Matrix IterationSolver::ArnoldiSolve(int m){
    if (m > A.rows)
        m = A.rows / 2;

    Matrix x = generateRandomMatrix(A.cols, 1, -1, 1);
    Matrix r = b - A * x;
    std::cout << "循环型阿诺尔迪算法开始！矩阵维度为: " << A.rows <<std::endl;
    std::cout << "初始残差为" << norm(r) << std::endl;
    
    int step = 1;
    while (norm(r) > epsilon){
        ArnoldResult ar_result = ArnoldiProcess(r, m);

        Matrix* v = new Matrix[m];
        std::copy_n(ar_result.v, m, v);
        Matrix H(m, m, 0);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                H(i, j) = ar_result.H(i, j);

        Matrix beta(m, 1, 0);
        beta(0, 0) = norm(r);

        // 利用列主元高斯消去法求解Hy=beta
        luResult result = DoolittleSolver::luDecompose(H, true);
        Matrix L = result.LU_mat.lower_tri();
        Matrix U = result.LU_mat.upper_tri();
        for (int i = 0; i < L.rows; ++i)
            L(i, i) = 1;  // 将L矩阵的对角线元素置为1
            Matrix y_ = DoolittleSolver::solveByTri(L, result.LU_P * beta, "low");  // 求解Ly_ = Pbeta
        Matrix y = DoolittleSolver::solveByTri(U, y_, "up");  // 求解Uy = y_

        Matrix V(A.rows, m, 0);
        for (int i = 0; i < V.rows; ++i)
            for (int j = 0; j < V.cols; ++j)
                V(i, j) = v[j](i, 0);
        Matrix z = V * y;

        x = x + z;
        r = b - A * x;

        std::cout << "第" << step << "步残差为" << norm(r) << std::endl;
        step += 1;
        delete[] v;  // 释放v的内存
    }
    return x;
}

Matrix IterationSolver::GMRESSolve(int m){
    if (m > A.rows)
        m = A.rows / 2;

    Matrix x = generateRandomMatrix(A.cols, 1, -1, 1);
    Matrix r = b - A * x;
    std::cout << "循环型广义极小残余(GMRES)算法开始！矩阵维度为: " << A.rows <<std::endl;
    std::cout << "初始残差为" << norm(r) << std::endl;

    int step = 1;
    while (norm(r) > epsilon){
        ArnoldResult ar_result = ArnoldiProcess(r, m);

        Matrix beta(m + 1, 1, 0);
        beta(0, 0) = norm(r);

        // 此处求解min||beta e1 - Hm y||的方式可以改进
        Matrix y = (ar_result.H.transpose() * ar_result.H).inv() * 
        ar_result.H.transpose() * beta;

        Matrix V(A.rows, m, 0);
        for (int i = 0; i < V.rows; ++i)
            for (int j = 0; j < V.cols; ++j)
                V(i, j) = ar_result.v[j](i, 0);
        Matrix z = V * y;

        x = x + z;
        r = b - A * x;

        std::cout << "第" << step << "步残差为" << norm(r) << std::endl;
        step += 1;
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

void Chapter03Demos::Demo3(int m){
    double delta = 1e-3;
    Matrix B(10, 10, 0);
    for (int i = 0; i < 10; ++i){
        if (i == 0)
            B(i, i + 1) = -1 + delta;
        else if (i == 9)
            B(i, i - 1) = -1 - delta;
        else {
            B(i, i - 1) = -1 - delta;
            B(i, i + 1) = -1 + delta;
        }
        B(i, i) = 4;
    }
    
    Matrix I(10, 10, 0);
    for (int i = 0; i < 10; ++i)
        I(i, i) = 1;
    Matrix Z(10, 10, 0);
    
    Matrix *A = nullptr;
    for (int i = 0; i < m; ++i){
        Matrix *R = nullptr;
        if (i == 0)
            R = new Matrix(B);
        else if (i == 1)
            R = new Matrix(-I);
        else
            R = new Matrix(Z);

        for (int j = 1; j < m; ++j){
            Matrix M;
            if (i == j)
                M = B;
            else if (j == i + 1 || j == i - 1)
                M = -I;
            else
                M = Z;
            Matrix temp = concat(*R, M, 1);
            delete R;
            R = new Matrix(temp);
        }

        if (i == 0)
            A = new Matrix(*R);  // 深拷贝，防止解引用A时将R删除掉
        else{
            Matrix temp = concat(*A, *R, 0);
            delete A;
            A = new Matrix(temp);
        }
    }

    Matrix f((*A).rows, 1, 1);
    Matrix d = *A * f;

    IterationSolver solver(*A, d, 10000, 1e-3);
    Matrix x = solver.ArnoldiSolve();
}

void Chapter03Demos::Demo31(int m){
    Matrix A1 = generatePositiveDefiniteMatrix(m);
    Matrix A2_ = generateRandomMatrix(m, m);
    Matrix A2 = A2_ * A2_; // 保证与A1数值大小类似
    Matrix b = generateRandomMatrix(m, 1, -10, 10);

    std::cout << "对基于伽辽金原理的大型数组迭代解法进行测试，A为随机对称正定矩阵，b为[-10,10]随机矩阵：" << std::endl;
    IterationSolver solver1(A1, b, 10000, 1e-3);
    // 共轭梯度法
    Matrix x0 = solver1.ConjugateGradientSolve();
    // 循环型阿诺尔迪算法
    Matrix x1 = solver1.ArnoldiSolve();
    // 循环型广义极小残余算法
    Matrix x2 = solver1.GMRESSolve();
}
