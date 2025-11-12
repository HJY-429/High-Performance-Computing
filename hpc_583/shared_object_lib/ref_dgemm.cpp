#include "ref_dgemm.hpp"
#include <stdexcept>

void dgemm(double a, const std::vector <std::vector <double >> &A, 
            const std::vector <std::vector <double >> &B, double b, 
            std::vector <std::vector <double >> &C) {

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
    

    // C = alpha*A*B + beta*C
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < p; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = a * sum + b * C[i][j];
        }
    }
}