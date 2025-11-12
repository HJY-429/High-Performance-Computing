#include "ref_dgemv.hpp"
#include <stdexcept>

void dgemv(double a, const std::vector<std::vector<double>>& A,
           const std::vector<double>& x, double b, std::vector<double>& y) {
    size_t m = A.size();
    if (m == 0) throw std::invalid_argument("Matrix A must have at least one row");
    
    size_t n = A[0].size();
    if (n != x.size()) throw std::invalid_argument("Matrix A columns must match vector x size");
    if (m != y.size()) throw std::invalid_argument("Matrix A rows must match vector y size");
    
    for (const auto& row : A) {
        if (row.size() != n) {
            throw std::invalid_argument("All rows in matrix A must have same number of columns");
        }
    }
    
    // alpha*Ax + beta*y
    for (size_t i = 0; i < m; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        y[i] = a * sum + b * y[i];
    }
}