#include "cavity.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "utils.hpp"

CavitySolver::CavitySolver(int Nx_, int Ny_, double Re_, double dt_, double U_)
    : Nx(Nx_), Ny(Ny_), Lx(1.0), Ly(1.0), Re(Re_), dt(dt_), U(U_)
{
    dx = Lx / Nx;
    dy = Ly / Ny;
    int Ntot = (Nx + 2) * (Ny + 2);  // 还包含幽灵网格

    u.assign(Ntot, 0.0);
    v.assign(Ntot, 0.0);
    p.assign(Ntot, 0.0);
    u_star.assign(Ntot, 0.0);
    v_star.assign(Ntot, 0.0);
    rhs.assign(Ntot, 0.0);
}

void CavitySolver::apply_velocity_bc() {
    // 顶壁（给定速度边界条件）
    for (int i = 0; i <= Nx + 1; ++i) {
        u[idx(i, Ny)] = this->U;
        v[idx(i, Ny)] = 0.0;
        u[idx(i, Ny + 1)] = this->U; // 镜像
        v[idx(i, Ny + 1)] = 0.0;
    }

    // 底壁
    for (int i = 0; i <= Nx + 1; ++i) {
        u[idx(i, 0)] = 0.0;
        v[idx(i, 0)] = 0.0;
        u[idx(i, 1)] = 0.0;
        v[idx(i, 1)] = 0.0;
    }

    // 左右壁
    for (int j = 0; j <= Ny + 1; ++j) {
        // 左
        u[idx(0, j)] = 0.0;
        v[idx(0, j)] = 0.0;
        u[idx(1, j)] = 0.0;
        v[idx(1, j)] = 0.0;
        // 右
        u[idx(Nx + 1, j)] = 0.0;
        v[idx(Nx + 1, j)] = 0.0;
        u[idx(Nx, j)] = 0.0;
        v[idx(Nx, j)] = 0.0;
    }
}

void CavitySolver::compute_predictor(){
    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            int id = idx(i, j);

            double uij = u[id];
            double vij = v[id];

            // 中心差分对流项
            double du_dx = (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * dx);
            double du_dy = (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * dy);
            double dv_dx = (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * dx);
            double dv_dy = (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * dy);

            double conv_u = uij * du_dx + vij * du_dy;
            double conv_v = uij * dv_dx + vij * dv_dy;

            // 粘性项 Laplacian
            double lap_u = (u[idx(i + 1, j)] - 2.0 * uij + u[idx(i - 1, j)]) / (dx * dx)
                         + (u[idx(i, j + 1)] - 2.0 * uij + u[idx(i, j - 1)]) / (dy * dy);
            double lap_v = (v[idx(i + 1, j)] - 2.0 * vij + v[idx(i - 1, j)]) / (dx * dx)
                         + (v[idx(i, j + 1)] - 2.0 * vij + v[idx(i, j - 1)]) / (dy * dy);

            u_star[id] = uij + dt * ( -conv_u + (1 / Re) * lap_u );
            v_star[id] = vij + dt * ( -conv_v + (1 / Re) * lap_v );
        }
    }
}

void CavitySolver::build_poisson_rhs(){
    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            int id = idx(i, j);
            double du_dx = (u_star[idx(i + 1, j)] - u_star[idx(i - 1, j)]) / (2.0 * dx);
            double dv_dy = (v_star[idx(i, j + 1)] - v_star[idx(i, j - 1)]) / (2.0 * dy);
            rhs[id] = (du_dx + dv_dy) / dt;
        }
    }
}

void CavitySolver::solve_poisson(int max_iter, double tol){
    std::vector<double> p_new = p;
    double idx2 = 1.0 / (dx * dx);
    double idy2 = 1.0 / (dy * dy);
    double coef = 2.0 * (idx2 + idy2);

    for (int it = 0; it < max_iter; ++it) {
        double max_err = 0.0;

        for (int j = 1; j <= Ny; ++j) {
            for (int i = 1; i <= Nx; ++i) {
                int id = idx(i, j);
                double pE = p[idx(i + 1, j)];
                double pW = p[idx(i - 1, j)];
                double pN = p[idx(i, j + 1)];
                double pS = p[idx(i, j - 1)];

                double rhs_loc = rhs[id];

                double p_old = p[id];
                double p_new_ij = ( (pE + pW) * idx2 + (pN + pS) * idy2 - rhs_loc ) / coef;
                p_new[id] = p_new_ij;
                max_err = std::max(max_err, std::fabs(p_new_ij - p_old));
            }
        }

        p.swap(p_new);

        // 压力 Neumann 边界
        for (int i = 0; i <= Nx + 1; ++i) {
            p[idx(i, 0)]      = p[idx(i, 1)];
            p[idx(i, Ny + 1)] = p[idx(i, Ny)];
        }
        for (int j = 0; j <= Ny + 1; ++j) {
            p[idx(0, j)]      = p[idx(1, j)];
            p[idx(Nx + 1, j)] = p[idx(Nx, j)];
        }
        // 固定一个点 p=0 去除常数自由度
        p[idx(0, 0)] = 0.0;

        if (max_err < tol) {
            std::cout << "泊松方程在 " << it << " 步收敛, err = "
                      << max_err << std::endl;
            break;
        }
    }
}

void CavitySolver::correct_velocity(){
    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            int id = idx(i, j);
            double dp_dx = (p[idx(i + 1, j)] - p[idx(i - 1, j)]) / (2.0 * dx);
            double dp_dy = (p[idx(i, j + 1)] - p[idx(i, j - 1)]) / (2.0 * dy);

            u[id] = u_star[id] - dt * dp_dx;
            v[id] = v_star[id] - dt * dp_dy;
        }
    }
}

double CavitySolver::velocity_residual() const{
    double res = 0.0;
    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            int id = idx(i, j);
            res += u[id] * u[id] + v[id] * v[id];
        }
    }
    return std::sqrt(res / (Nx * Ny));
}

void CavitySolver::write_output(const std::string& filename) const{
    std::ofstream ofs(filename);
    ofs << "x,y,u,v,p\n";
    for (int j = 1; j <= Ny; ++j) {
        for (int i = 1; i <= Nx; ++i) {
            int id = idx(i, j);
            double x = (i - 0.5) * dx;
            double y = (j - 0.5) * dy;
            ofs << x << "," << y << ","
                << u[id] << "," << v[id] << ","
                << p[id] << "\n";
        }
    }
    std::cout << "数据写入: " << filename << std::endl;
}

void CavitySolver::run(double t_end, int output_interval_steps, int poisson_iter, double poisson_tol){
    int    step = 0;
    double t    = 0.0;

    while (t < t_end) {
        // --- 投影法 时间步 ---
        apply_velocity_bc();       // 旧时间层/边界
        compute_predictor();       // Step 1: u*, v*
        apply_velocity_bc();       // 对 u* 再施边界
        build_poisson_rhs();       // Step 2: div(u*)/dt
        solve_poisson(poisson_iter, poisson_tol); // Step 2: 解 p^{n+1}
        correct_velocity();        // Step 3: u^{n+1}, v^{n+1}
        apply_velocity_bc();       // 保证边界速度

        // 打印信息方便看收敛/进度
        if (step % 100 == 0) {
            double res = velocity_residual();
            std::cout << "step = " << step
                      << ", t = " << t
                      << ", vel L2 ~ " << res << std::endl;
        }

        // --- 每隔 output_interval_steps 步写一次文件 ---
        if (output_interval_steps > 0 && (step % output_interval_steps == 0)) {
            std::ostringstream name;
            // 例如：cavity_step000000_t0.000000.csv
            name << "cavity_data/cavity_step"
                 << std::setw(6) << std::setfill('0') << step
                 << "_t" << std::fixed << std::setprecision(6) << t
                 << ".csv";
            write_output(name.str());
        }

        // 更新时间 & 步数
        ++step;
        t += dt;

        // 防止由于 dt 不整除 t_end 导致无限循环
        if (t + 0.5*dt > t_end) {
            t = t_end;
        }
    }

    // 最后一帧再写一次（确保有 t_end 的结果）
    std::ostringstream name;
    name << "cavity_data/cavity_final_t" << std::fixed << std::setprecision(6) << t_end << ".csv";
    write_output(name.str());
}
