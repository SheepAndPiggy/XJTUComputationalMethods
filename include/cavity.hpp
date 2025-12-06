#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "utils.hpp"

struct CavitySolver{
    int Nx, Ny;  // x和y方向的网格数量
    double Lx, Ly;  // x和y方向的特征长度
    double dt;  // 时间步长
    double dx, dy;  // 网格步长，为Lx/Nx,Ly/Ny
    double U;  // 顶部速度
    double Re;  // 雷诺数 ULx/nu

    std::vector<double> u, v, p;  // 当前时间层的速度和压力
    std::vector<double> u_star, v_star;  // 投影法中间预测速度
    std::vector<double> rhs;  // 压力泊松方程右侧项

    // 构造函数
    CavitySolver(int Nx_, int Ny_, double nu_, double dt_, double U_); 

    // 将二维的索引转换成一维的索引，方便获取数据
    inline int idx(int i, int j) const {
        // i = 0..Nx+1, j = 0..Ny+1
        return i + (Nx + 2) * j;
    }

    // 设置速度边界条件
    void apply_velocity_bc();

    // 计算中间速度（预测速度）
    void compute_predictor();

    // 构造压力泊松方程的右端
    void build_poisson_rhs();

    // 解泊松方程（雅可比迭代法）
    void solve_poisson(int max_iter, double tol);

    // 将中间速度投影，使其满足连续方程
    void correct_velocity();

    // 速度场L2范数（判断是否收敛）
    double velocity_residual() const;

    // 将结果写入csv文件
    void write_output(const std::string& filename) const;

    // 运行函数
    void run(double t_end,
         int    output_interval_steps,
         int    poisson_iter = 5000,
         double poisson_tol  = 1e-6);
};