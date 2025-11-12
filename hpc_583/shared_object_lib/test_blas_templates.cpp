#include <iostream>
#include <vector>
#include "ref_axpyt.hpp"
#include "ref_gemvt.hpp"
#include "ref_gemmt.hpp"

template <typename T>
void test_axpy() {
    std::vector<T> x = {1, 2, 3};
    std::vector<T> y = {4, 5, 6};
    T alpha = 2;
    
    axpy(alpha, x, y);
    
    std::cout << "axpy result: ";
    for (const auto& val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void test_gemv() {
    std::vector<std::vector<T>> A = {{1, 3}, {3, 4}, {5, 6}};
    std::vector<T> x = {2, 3};
    std::vector<T> y = {1, 1, 1};
    T alpha = 2, beta = 0.5;
    
    gemv(alpha, A, x, beta, y);
    
    std::cout << "gemv result: ";
    for (const auto& val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void test_gemm() {
    std::vector<std::vector<T>> A = {{1, 2}, {3, 4}};
    std::vector<std::vector<T>> B = {{5, 6}, {7, 8}};
    std::vector<std::vector<T>> C = {{1, 0}, {0, 1}};
    T alpha = 2, beta = 0.5;
    
    gemm(alpha, A, B, beta, C);
    
    std::cout << "gemm result:\n";
    for (const auto& row : C) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    std::cout << "Testing with double:\n";
    test_axpy<double>();
    test_gemv<double>();
    test_gemm<double>();
    
    std::cout << "\nTesting with float:\n";
    test_axpy<float>();
    test_gemv<float>();
    test_gemm<float>();

    std::cout << "\nTesting with int:\n";
    test_axpy<int>();
    test_gemv<int>();
    test_gemm<int>();
    
    return 0;
}