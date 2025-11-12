#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cblas.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>

void init_matrix(double* mat, int m, int n) {
    for (int i = 0; i < m * n; i++) {
        mat[i] = (double)rand() / RAND_MAX;
    }
}

void openBLAS_perform(int min_n, int max_n, int ntrials, double alpha, double beta, const std::string& filename){
    // Create and open a CSV file
    std::ofstream outfile(filename);
    outfile << "n,Avg_Time(s),FLOPs\n";

    for (int n = min_n; n <= max_n; n *= 2){
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
        
        avg_time = elapsed_time / ntrials;
        double flops = (2.0 * n * n * n + 2.0 * n * n) / avg_time;
        
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        outfile << n << "," << avg_time << "," << flops << "\n";

        delete[] A;
        delete[] B;
        delete[] C;
    }
    outfile.close();
}

void CUBLAS_perform(int min_n, int max_n, int ntrials, double alpha, double beta, const std::string& filename){
    // Create and open a CSV file
    std::ofstream outfile(filename);
    outfile << "n,Avg_Time(s),FLOPs\n";

    cublasHandle_t handle;
    if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "Failed to initialize CUBLAS" << std::endl;
        return;
    }

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    for (int n = min_n; n <= max_n; n *= 2){
        double* hA = new double[n * n];
        double* hB = new double[n * n];
        double* hC = new double[n * n];
        double *dA, *dB, *dC;

        // initialize
        double elapsed_time = 0.0;
        double avg_time;

        init_matrix(hA, n, n);
        init_matrix(hB, n, n);
        init_matrix(hC, n, n);

        cudaMalloc(&dA, n * n * sizeof(double));
        cudaMalloc(&dB, n * n * sizeof(double));
        cudaMalloc(&dC, n * n * sizeof(double));

        for (int t = 0; t < ntrials; t++) {
            cudaMemcpy(dA, hA, n * n * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(dB, hB, n * n * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(dC, hC, n * n * sizeof(double), cudaMemcpyHostToDevice);

            cudaEventRecord(start);
            
            cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        n, n, n, &alpha, dA, n, dB, n, &beta, dC, n);
            
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);
            elapsed_time += milliseconds;
        }

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);
        
        avg_time = (elapsed_time / ntrials) * 1.e-3;
        double flops = (2.0 * n * n * n + 2.0 * n * n) / avg_time;
        
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        outfile << n << "," << avg_time << "," << flops << "\n";
    }

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cublasDestroy(handle);
    outfile.close();
}

int main(){
    const int ntrials = 3;
    const int min_n = 2;
    const int max_n = 16384;
    const double alpha = 2; 
    const double beta = 0.5;

    openBLAS_perform(min_n, max_n, ntrials, alpha, beta, "openBLAS_performance.csv");
    CUBLAS_perform(min_n, max_n, ntrials, alpha, beta, "CUBLAS_performance.csv");


    return 0;
}
