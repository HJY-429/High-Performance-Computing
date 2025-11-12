#include "ref_daxpy.hpp"
#include <stdexcept>

void daxpy(double a, const std::vector <double > &x, std::vector <double > &y){
    if (x.size() != y.size()){
        throw std::invalid_argument("Vectors sizes must be equal!");
    }
    if (x.empty() || y.empty()) {
        throw std::invalid_argument("Vectors should not be empty");
    }

    for (size_t i=0; i < x.size(); ++i){
        y[i] += a * x[i];
    }
}