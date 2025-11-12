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
#include <lapacke.h>

double compute_residual(std::complex<double>* A, std::complex<double>* b, std::complex<double>* z, int n) {
    std::vector<std::complex<double>> Ax(n, 0);
    std::complex<double> alpha(1.0, 0.0);
    std::complex<double> beta(0.0, 0.0);
    
    cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, &alpha, 
                A, n, z, 1, &beta, Ax.data(), 1);
    
    for (int i = 0; i < n; i++) {
        Ax[i] = b[i] - Ax[i];
    }
    
    return cblas_dznrm2(n, Ax.data(), 1);
}

double compute_infinity_norm(std::complex<double>* A, int n) {
    double max = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += std::abs(A[i + j*n]);
        }
        if (row_sum > max) max = row_sum;
    }
    return max;
}

int main() {
    const double epsilon = std::numeric_limits<double>::epsilon();
    std::vector<int> sizes = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
    
    // Open CSV file for writing
    std::ofstream ofs("solver_data.csv");
    ofs << "matrixSize,residual,log10residual,normalizedError,log10normalizedError\n";
    
    for (int n : sizes) {
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
        
        // Compute metrics
        double residual = compute_residual(A_orig.data(), b_orig.data(), z, n);
        double A_inf_norm = compute_infinity_norm(A_orig.data(), n);
        double z_norm = cblas_dznrm2(n, z, 1);
        double normal_error = residual / (A_inf_norm * z_norm * epsilon);
        
        // Write to CSV
        ofs << n << "," << residual << "," << log10(residual) << "," 
        << normal_error << "," << log10(normal_error) << "\n";
        
        free(A);
        free(b);
        free(z);
    }
    
    ofs.close();
    std::cout << "Data saved as 'solver_data.csv'" << std::endl;
    
    return 0;
}