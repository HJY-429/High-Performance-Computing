#include "ref_gemvt.hpp"
#include <stdexcept>

template <typename T>
void gemv(T a, const std::vector<std::vector<T>>& A, const std::vector<T>& x, T b, std::vector<T>& y) {
  // Check dimensions
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

  for (size_t i = 0; i < m; ++i) {
    T temp = T();
    for (size_t j = 0; j < n; ++j) {
      temp += A[i][j] * x[j];
    }
    y[i] = a * temp + b * y[i];
  }
}

template void gemv(int a, const std::vector<std::vector<int>>& A, 
  const std::vector<int>& x, int b, std::vector<int>& y);
template void gemv(float a, const std::vector<std::vector<float>>& A, 
  const std::vector<float>& x, float b, std::vector<float>& y);
template void gemv(double a, const std::vector<std::vector<double>>& A, 
  const std::vector<double>& x, double b, std::vector<double>& y);
