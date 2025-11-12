#ifndef REF_GEMMT_HPP
#define REF_GEMMT_HPP

#include <vector>

template <typename T>
void gemm(T alpha, const std::vector<std::vector<T>>& A,
          const std::vector<std::vector<T>>& B, T beta,
          std::vector<std::vector<T>>& C);

#endif