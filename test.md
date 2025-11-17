```c++
Matrix IterationSolver::GradientSolve(){
    Matrix x = generateRandomMatrix(A.cols, 1, -1, 1);
    
    Matrix gradient = A * x - b;
    Matrix r = b - A * x;

    for (int i = 0; i < max_iter; ++i){
        double final_error = this->residualNormError(x) / norm(x);
        if (final_error < epsilon){
            std::cout << "维度为" << A.cols << "的梯度下降法过程共用" << i + 1 << "步, "
             << "残差||Ax-b||/||x||为" << final_error << std::endl;
            break;
        }
        double alpha = (r.transpose() * r)(0, 0) / (r.transpose() * A * r)(0, 0);
        x = x - alpha * gradient;
        gradient = A * x - b;
        r = b - A * x;
    }
    return x;
}
```