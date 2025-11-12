#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include "ref_dgemm.hpp"

void init_matrix(std::vector<std::vector<double>>& A, int m, int n) {
    A.resize(m, std::vector<double>(n));
    for (auto& row : A) {
        for (auto& elem : row) {
            elem = 1.0;
        }
    }
}

void dgemm_perform(){
    const int ntrials = 3;
    const int min_n = 2;
    const int max_n = 512;
    
    // Create and open a CSV file
    std::ofstream outfile("dgemm_performance.csv");
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n++){
        std::vector<std::vector<double>> A, B, C;

        // initialize
        double alpha = 2; double beta = 0.5;
        long double elapsed_time = 0.L;
        long double avg_time;

        init_matrix(A, n, n);
        init_matrix(B, n, n);
        init_matrix(C, n, n);

        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            dgemm(alpha, A, B, beta, C);
            
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

void test_dgemm() {
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}};
    std::vector<std::vector<double>> B = {{5, 6}, {7, 8}};
    std::vector<std::vector<double>> C = {{1, 0}, {0, 1}};
    double alpha = 2, beta = 0.5;
    
    dgemm(alpha, A, B, beta, C);
    
    std::cout << "gemm result:\n";
    for (const auto& row : C) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void test_err_dgemm() {
    std::vector<std::vector<double>> Ae = {{1, 2}, {3, 4}};
    std::vector<std::vector<double>> Be = {{5, 6, 7}, {7, 8, 9}};
    std::vector<std::vector<double>> Ce = {{1, 0}, {0, 1}};
    double alpha = 2, beta = 0.5;
    
    dgemm(alpha, Ae, Be, beta, Ce);
    
    std::cout << "gemm result:\n";
    for (const auto& row : Ce) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

int main(){
    dgemm_perform();
    test_dgemm();
    // test_err_dgemm();
        
}
