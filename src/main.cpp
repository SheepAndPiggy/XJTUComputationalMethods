#include <iostream>
#include "utils.hpp"
#include "chapter02.hpp"
#include "chapter03.hpp"
#include "chapter04.hpp"
#include "chapter05.hpp"
#include "chapter07.hpp"
#include "cavity.hpp"
#include <windows.h>

void cavity_demo();

int main(){
    SetConsoleOutputCP(CP_UTF8);  // 强制控制台使用utf-8编码
    cavity_demo();
    return 0;
}

void cavity_demo() {
    int    Nx = 64;
    int    Ny = 64;
    double Re = 1000;
    double dt = 0.001;
    double U = 1.0;

    CavitySolver solver(Nx, Ny, Re, dt, U);

    double t_end = 5.0;          // 仿真总时长（物理时间）
    int    output_interval = 500; // 每 500 步输出一次

    solver.run(t_end, output_interval);
}

