#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include "ref_dgemv.hpp"

void init_matrix(std::vector<std::vector<double>>& A, int m, int n) {
    A.resize(m, std::vector<double>(n));
    for (auto& row : A) {
        for (auto& elem : row) {
            elem = 1.0;
        }
    }
}

void init_vec(std::vector<double>& v, int size) {
    v.resize(size);
    for (auto& elem : v) {
        elem = 1.0;
    }
}

void dgemv_perform(){
    const int ntrials = 3;
    const int min_n = 2;
    const int max_n = 512;
    
    // Create and open a CSV file
    std::ofstream outfile("dgemv_performance.csv");
    outfile << "n,Avg_Time(s),FLOPs\n";
    

    for (int n = min_n; n <= max_n; n++){
        std::vector<std::vector<double>> A;
        std::vector<double> x, y;

        // initialize
        double alpha = 2; double beta = 0.5;
        long double elapsed_time = 0.L;
        long double avg_time;

        init_matrix(A, n, n);
        init_vec(x, n);
        init_vec(y, n);

        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            dgemv(alpha, A, x, beta, y);
            
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

void test_dgemv(){
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}, {5, 6}};
    std::vector<double> x = {2, 3};
    std::vector<double> y = {1, 1, 1};
    double alpha = 2, beta = 0.5;
    
    dgemv(alpha, A, x, beta, y);
    
    std::cout << "gemv result: ";
    for (const auto& val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void test_err_dgemv(){
    std::vector<std::vector<double>> Ae1 = {{1, 2}, {3, 4}, {5, 6, 7}};
    std::vector<double> xe = {2, 3};
    std::vector<double> ye = {1, 1, 1};
    double alpha = 2, beta = 0.5;
    
    dgemv(alpha, Ae1, xe, beta, ye);
    
    std::cout << "gemv result: ";
    for (const auto& val : ye) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

int main(){
    dgemv_perform();
    test_dgemv();
    // test_err_dgemv();
}
