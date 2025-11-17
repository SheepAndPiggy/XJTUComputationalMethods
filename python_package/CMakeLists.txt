cmake_minimum_required(VERSION 3.12)
project(CM_pkg)

# 设置 pybind11 的目录路径
set(pybind11_DIR "C:/Users/dan/AppData/Local/Programs/Python/Python313/Lib/site-packages/pybind11/share/cmake/pybind11")

# 寻找 pybind11 包
set(PYBIND11_FINDPYTHON ON)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# 设置静态链接
set(BUILD_SHARED_LIBS OFF)

# 添加 C++ 模块
pybind11_add_module(CM_pkg CM_pkg.cpp) # 模块名称 'CM_pkg'

# 设置输出目录
set_target_properties(CM_pkg PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# 针对 Windows 平台的额外设置
if (WIN32)
    target_compile_definitions(CM_pkg PRIVATE NOMINMAX)
    target_link_libraries(CM_pkg PRIVATE ${PYTHON_LIBRARIES})

    # 针对 MinGW 添加静态链接选项
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_link_options(CM_pkg PRIVATE -static-libgcc -static-libstdc++)
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(CM_pkg PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
    target_compile_options(CM_pkg PRIVATE /W4 /permissive-)
endif()

message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Using Python include directories: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Output library will be placed in: ${CMAKE_CURRENT_SOURCE_DIR}")