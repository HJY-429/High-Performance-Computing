#include "refBLAS.hpp"
#include <iostream>
#include <vector>

void test_daxpy() {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {4.0, 5.0, 6.0};
    daxpy(2.0, x, y);
    
    std::cout << "daxpy result: ";
    for (const auto& val : y) std::cout << val << " ";
    std::cout << std::endl;
}

void test_dgemv() {
    std::vector<std::vector<double>> A = {{1.0, 2.0}, {3.0, 4.0}};
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> y = {0.0, 0.0};
    dgemv(1.0, A, x, 0.0, y);
    
    std::cout << "dgemv result: ";
    for (const auto& val : y) std::cout << val << " ";
    std::cout << std::endl;
}

void test_dgemm() {
    std::vector<std::vector<double>> A = {{1.0, 2.0}, {3.0, 4.0}};
    std::vector<std::vector<double>> B = {{5.0, 6.0}, {7.0, 8.0}};
    std::vector<std::vector<double>> C = {{1.0, 0.0}, {1.0, 0.0}};
    dgemm(2.0, A, B, 0.5, C);
    
    std::cout << "dgemm result:\n";
    for (const auto& row : C) {
        for (const auto& val : row) std::cout << val << " ";
        std::cout << std::endl;
    }
}

void test_templates() {
    std::vector<float> xf = {1.0f, 2.0f};
    std::vector<float> yf = {3.0f, 4.0f};
    axpy(2.0f, xf, yf);
    
    std::cout << "float axpy result: ";
    for (const auto& val : yf) std::cout << val << " ";
    std::cout << std::endl;
}

int main() {
    test_daxpy();
    test_dgemv();
    test_dgemm();
    test_templates();
    return 0;
}
