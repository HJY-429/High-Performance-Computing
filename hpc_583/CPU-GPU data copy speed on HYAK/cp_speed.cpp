#include <iostream>
#include <cuda.h>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cuda_runtime.h>

double measure_bw(size_t size, bool h2d){

    char* dataH = new char [size];
    char* dataD;
    cudaMalloc(&dataD, size);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    
    for (size_t i = 0; i < size; i++){
        dataH[i] = i % 256; 
    }
    cudaEventRecord(start);
    if (h2d){
        cudaMemcpy(dataD, dataH, size, cudaMemcpyHostToDevice);
    } else{
        cudaMemcpy(dataH, dataD, size, cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    cudaFree(dataD);
    delete[] dataH;
    return (milliseconds * 1.e-3);
}

double count_bw(double t, size_t bytes){
    if (t<=0) return 0;
    double GB = static_cast<double>(bytes) / (1024 * 1024 * 1024);
    return (GB / t);
}

int main(){
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0){
        std::cerr << "No CUDA device found" << std::endl;
        return 1;
    }

    cudaSetDevice(0);
    std::vector<size_t> sizes;
    for (size_t size = 1; size <= 2ULL * 1024 * 1024 * 1024; size *=2){
        sizes.push_back(size);
    }

    std::ofstream csv_file("copy_bw.csv");
    csv_file << "size_bytes,size_mb,h2d_time,h2d_bw,d2h_time,d2h_bw\n";

    std::cout << std::setw(12) << "Size (B)" << std::setw(12) << "Size (MB)" << std::setw(12) << "H2D Time" 
    << std::setw(12) << "H2D BW" << std::setw(12) << "D2H Time" << std::setw(12) << "D2H BW" << std::endl;

    for (size_t size : sizes){
        double size_mb = static_cast<double>(size) / (1024 * 1024);

        // H2D performance
        double h2d_t = measure_bw(size, true);
        double h2d_bw = count_bw(h2d_t, size);
        // D2H performance
        double d2h_t = measure_bw(size, false);
        double d2h_bw = count_bw(d2h_t, size);

        std::cout << std::setw(10) << size << "B" << std::setw(10) << size_mb << "MB" 
        << std::setw(10) << h2d_t << "s" << std::setw(10) << h2d_bw << "GB/s" 
        << std::setw(10) << d2h_t << "s" << std::setw(10) << d2h_bw << "GB/s" << std::endl;

        csv_file << size << "," 
                 << size_mb << "," 
                 << h2d_t << "," << h2d_bw << "," 
                 << d2h_t << "," << d2h_bw << "\n";
    }

    csv_file.close();
    std::cout << "Results saved as 'copy_bw.csv'" << std::endl;

    return 0;
}