#include "ref_gemmt.hpp"
#include <stdexcept>

template <typename T>
void gemm(T alpha, const std::vector<std::vector<T>>& A,
          const std::vector<std::vector<T>>& B, T beta,
          std::vector<std::vector<T>>& C) {
    // Check if matrices are empty
    if (A.empty() || A[0].empty()) throw std::invalid_argument("Matrix A should not be empty");
    if (B.empty() || B[0].empty()) throw std::invalid_argument("Matrix B should not be empty");

    const size_t m = A.size();    // Rows of A
    const size_t p = A[0].size(); // Columns of A / Rows of B
    const size_t n = B[0].size(); // Columns of B

    // Verify A and B are properly rectangulars
    for (const auto& row : A) {
        if (row.size() != p) {
            throw std::invalid_argument("Matrix A has inconsistent row sizes");
        }
    }
    for (const auto& row : B) {
        if (row.size() != n) {
            throw std::invalid_argument("Matrix B has inconsistent row sizes");
        }
    }
    // Verify A's columns match B's rows
    if (B.size() != p) {
        throw std::invalid_argument(
            "Matrix multiplication dimension mismatch: A columns (" + 
            std::to_string(p) + ") != B rows (" + std::to_string(B.size()) + ")"
        );
    }
    
    // Verify C has correct dimensions
    if (C.size() != m) {
        throw std::invalid_argument(
            "Output matrix C has wrong number of rows: expected " + 
            std::to_string(m) + ", got " + std::to_string(C.size())
        );
    }
    for (const auto& row : C) {
        if (row.size() != n) {
            throw std::invalid_argument(
                "Output matrix C has wrong number of columns: expected " + 
                std::to_string(n) + ", got " + std::to_string(row.size())
            );
        }
    }
    
    // alpha*AB + beta*C
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            T temp = T(); //T()
            for (size_t k = 0; k < p; ++k) {
                temp += A[i][k] * B[k][j];
            }
            C[i][j] = alpha * temp + beta * C[i][j];
        }
    }
}

template void gemm(int alpha, const std::vector<std::vector<int>>& A,
    const std::vector<std::vector<int>>& B, int beta,
    std::vector<std::vector<int>>& C);
template void gemm(float alpha, const std::vector<std::vector<float>>& A,
    const std::vector<std::vector<float>>& B, float beta,
    std::vector<std::vector<float>>& C);
template void gemm(double alpha, const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B, double beta,
    std::vector<std::vector<double>>& C);
