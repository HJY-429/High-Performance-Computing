#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <cblas.h>

// void init_matrix(std::vector<double>& A, int rows, int cols, double value=1.0) {
//     A.resize(rows * cols);
//     for (auto& elem : A) {
//             elem = value;
//         }
// }

void init_matrix(double* mat, int m, int n) {
    for (int i = 0; i < m * n; i++) {
        mat[i] = (double)rand() / RAND_MAX;
    }
}

// void init_vector(std::vector<double>& x, int n, double value=1.0) {
//     x.resize(n);
//     for (auto& elem : x) {
//             elem = value;
//         }
// }

void init_vector(double* vec, int n) {
    for (int i = 0; i < n; i++) {
        vec[i] = (double)rand() / RAND_MAX;
    }
}

void daxpy_perform(int min_n, int max_n, int ntrials, double alpha, double beta, const std::string& filename){
    // Create and open a CSV file
    std::ofstream outfile(filename);
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n*=2){
        // std::vector<double> x, y;
        double* x = new double[n];
        double* y = new double[n];

        // initialize
        long double elapsed_time = 0.L;
        long double avg_time;

        init_vector(x, n);
        init_vector(y, n);

        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            cblas_daxpy(n, alpha, x, 1, y, 1);
            
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        
        avg_time = elapsed_time / static_cast<long double>(ntrials);
        double flops = (2.0 * n) / avg_time;
        
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        outfile << n << "," << avg_time << "," << flops << "\n";
    }
    outfile.close();
}

void dgemv_perform(int min_n, int max_n, int ntrials, double alpha, double beta, const std::string& filename){
    // Create and open a CSV file
    std::ofstream outfile(filename);
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n*=2){
        double* A = new double[n * n];
        double* x = new double[n];
        double* y = new double[n];

        // initialize
        long double elapsed_time = 0.L;
        long double avg_time;

        init_matrix(A, n, n);
        init_vector(x, n);
        init_vector(y, n);

        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, beta, y, 1);
            
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        
        avg_time = elapsed_time / static_cast<long double>(ntrials);
        double flops = (2.0 * n * n + 2.0 * n) / avg_time;
        
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        outfile << n << "," << avg_time << "," << flops << "\n";
    }
    outfile.close();
}

void dgemm_perform(int min_n, int max_n, int ntrials, double alpha, double beta, const std::string& filename){
    // Create and open a CSV file
    std::ofstream outfile(filename);
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n*=2){
        double* A = new double[n * n];
        double* B = new double[n * n];
        double* C = new double[n * n];

        // initialize
        long double elapsed_time = 0.L;
        long double avg_time;

        init_matrix(A, n, n);
        init_matrix(B, n, n);
        init_matrix(C, n, n);

        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        n, n, n, alpha, A, n, B, n, beta, C, n);
            
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        
        avg_time = elapsed_time / static_cast<long double>(ntrials);
        double flops = (2.0 * n * n * n + 2.0 * n * n) / avg_time;
        
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        outfile << n << "," << avg_time << "," << flops << "\n";
    }
    outfile.close();
}

int main(){
    const int ntrials = 3;
    const int min_n = 2;
    const int max_n = 4096;
    double alpha = 2.0; 
    double beta = 0.5;

    daxpy_perform(min_n, max_n, ntrials, alpha, beta, "daxpy_performance_t.csv");
    dgemv_perform(min_n, max_n, ntrials, alpha, beta, "dgemv_performance_t.csv");
    dgemm_perform(min_n, max_n, ntrials, alpha, beta, "dgemm_performance_t.csv");

    return 0;
}
