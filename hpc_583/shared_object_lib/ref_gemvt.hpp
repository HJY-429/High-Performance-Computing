#ifndef REF_GEMVT_HPP
#define REF_GEMVT_HPP

#include <vector>

template <typename T>
void gemv(T a, const std::vector<std::vector<T>>& A,
          const std::vector<T>& x, T b, std::vector<T>& y);

#endif 