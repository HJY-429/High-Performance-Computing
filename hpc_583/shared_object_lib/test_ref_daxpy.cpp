#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include "ref_daxpy.hpp"

void test_daxpy(){
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5, 6};
    double alpha = 2;
    
    daxpy(alpha, x, y);
    
    std::cout << "axpy result: ";
    for (const auto& val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void test_err_daxpy(){
    std::vector<double> x = {1, 2, 3};
    std::vector<double> y = {4, 5};
    double alpha = 2;
    
    daxpy(alpha, x, y);
    
    std::cout << "axpy result: ";
    for (const auto& val : y) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void daxpy_perform(){
    const int ntrials = 3;
    const int min_n = 2;
    const int max_n = 512;

    // Create and open a CSV file
    std::ofstream outfile("daxpy_performance.csv");
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n++){
        std::vector<double> x(n), y(n);
        double alpha = 2;
        long double elapsed_time = 0.L;
        long double avg_time;

        for (int i = 0; i < n; i++) {
            x[i] = 1.0; y[i] = 1.0;
        }
        for (int t = 0; t < ntrials; t++) {
            auto start = std::chrono::high_resolution_clock::now();
            
            daxpy(alpha, x, y);
            
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
int main(){

    daxpy_perform();
    test_daxpy();
    // test_err_daxpy();
        
}
