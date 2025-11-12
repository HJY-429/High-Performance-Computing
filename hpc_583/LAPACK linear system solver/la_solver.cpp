#include <iostream>
#include <complex>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#include <limits>
#include <fstream>
#include <cblas.h>
#include <complex.h>
#include <lapacke.h>

void cgemv(int n, std::complex<double>* A, std::complex<double>* x, std::complex<double>* y){
    for (int i = 0; i < n; i++){
        y[i] = 0;
        for (int j = 0; j < n; j++){
            y[i] += A[i + j * n] * x[j];
        }
    }
}

// compute complex 2-norm
double complex_norm(std::complex<double>* x, int n){
    double norm = 0.0;
    for (int i = 0; i < n; i++){
        norm += std::norm(x[i]);
    }
    return std::sqrt(norm);
}

// compute || b - Az ||
double comp_residual(std::complex<double>* A, std::complex<double>* b, std::complex<double>* z, int n){
    std::vector<std::complex<double>> Az(n);
    cgemv(n, A, z, Az.data());
    for (int i = 0; i < n; i++){
        Az[i] = b[i] - Az[i];
    }

    return complex_norm(Az.data(), n);
}

double comp_inf_norm(std::complex<double>* A, int n){
    double max = 0.0;
    for(int i = 0; i < n; i++){
        double row_sum = 0.0;
        for (int j = 0; j < n; j++){
            row_sum += std::abs(A[i + j * n]);
        }
        if (row_sum > max)  max = row_sum;
    }
    return max;
}

int main() {
    const double eps = std::numeric_limits<double>::epsilon();
    std::vector<int> sizes;
    for (int size = 16; size <= 8192; size *=2){
        sizes.push_back(size);
    }

    std::ofstream ofs("la_solver.csv");
    ofs << "matrixSize,residual,log10residual,normalizedError,log10normalizedError\n";
    for (int n : sizes){
        std::complex<double>* A = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n * n); 
        std::complex<double>* b = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n);
        std::complex<double>* z = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n);
        std::vector<std::complex<double>> A_orig(n * n);
        std::vector<std::complex<double>> b_orig(n);
        std::vector<int> ipiv(n);

        srand (0);
        int k =0;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n ; i++) {
                A[k] = 0.5 - (double)rand() / (double)RAND_MAX + std::complex <double >(0, 1)
                * (0.5 - (double)rand() / (double)RAND_MAX); 
                if (i == j)     A[k] *= static_cast<double>(n);
                A_orig[k] = A[k];
                k++;
            } 
        }
        srand (1);
        for (int i = 0; i < n; i++) {
            b[i] = 0.5 - (double)rand() / (double)RAND_MAX + std::complex <double >(0, 1)
            * (0.5 - (double)rand() / (double)RAND_MAX);
            z[i] = b[i];
            b_orig[i] = b[i];
        }

        int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, 1, reinterpret_cast<lapack_complex_double*>(A), n, 
                    ipiv.data(), reinterpret_cast<lapack_complex_double*>(z), n);
        if (info != 0){
            std::cerr << "Solution failed for n = " << n << std::endl;
            continue;
        }

        double residual = comp_residual(A_orig.data(), b_orig.data(), z, n);
        double A_infNorm = comp_inf_norm(A_orig.data(), n);
        double z_norm = complex_norm(z, n);

        double normal_error = residual / (A_infNorm * z_norm * eps);

        ofs << n << "," << residual << "," << log10(residual) << "," 
        << normal_error << "," << log10(normal_error) << "\n";

        free(A);
        free(b);
        free(z);
    }

    ofs.close();
    std::cout << "Data saved as 'la_solver.csv'" << std::endl;

    return 0;
}