#include <iostream>
#include "utils.hpp"
#include "chapter02.hpp"
#include "chapter03.hpp"
#include <windows.h>

int main(){
    SetConsoleOutputCP(CP_UTF8);  // 强制控制台使用utf-8编码
    Chapter03Demos::Demo31(1000);
    return 0;
}
