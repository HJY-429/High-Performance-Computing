#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

// write
void matrix_wr(const std::string& filename, const std::vector<double>& matrix, int n){
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs){
        std::cerr << "Failed to open file for writing.\n";
    }
    // column major matrix
    for (int col = 0; col < n; ++col){
        for (int row = 0; row < n; ++row){
            ofs.write(reinterpret_cast<const char*>(&matrix[row + col * n]), sizeof(double));
        }
    }  
}

// read
void matrix_rd(const std::string& filename, std::vector<double>& matrix, int n){
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs){
        std::cerr << "Failed to open file for reading.\n";
    }
    // column major matrix
    for (int col = 0; col < n; ++col){
        for (int row = 0; row < n; ++row){
            ifs.read(reinterpret_cast<char *>(&matrix[row + col * n]), sizeof(double));
        }
    }  
}

double measure(int n, bool write){
    std::vector<double> matrix(n * n);

    for (int i = 0; i < n * n; ++i) {
        matrix[i] = static_cast<double>(i);
    }
    
    std::string filename = "matrix_" + std::to_string(n) + ".bin";
    auto start = std::chrono::high_resolution_clock::now();
    if (write){
        matrix_wr(filename, matrix, n);
    } else{
        matrix_rd(filename, matrix, n);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast <std::chrono::nanoseconds>(stop - start);

    return (duration.count() * 1.e-9);
}

double count_bw(double t, size_t bytes){
    if (t<=0) return 0;
    double GB = static_cast<double>(bytes) / (1024 * 1024 * 1024);
    return (GB / t);
}


void clean_files(const std::vector<std::string>& filenames) {
    std::cout << "\nCleaning up temporary files...\n";
    for (const auto& filename : filenames) {
        try {
            if (fs::remove(filename)) {
                std::cout << "Deleted: " << filename << "\n";
            }
        } catch (const fs::filesystem_error& err) {
            std::cerr << "Error deleting " << filename << ": " << err.what() << "\n";
        }
    }
}

int main(){
    std::vector<int> dim;
    for (int i = 5; i <= 14; ++i){
        dim.push_back(1 << i);
    }

    std::vector<std::string> binaryFiles;
    std::ofstream csv_file("matrix_bw.csv");
    csv_file << "dimension,matrix_size_mb,write_time,write_bw,read_time,read_bw\n";

    for (int n : dim){
        size_t matrix_bytes = n * n * sizeof(double);
        double matrix_mb = static_cast<double>(matrix_bytes) / (1024 * 1024);

        std::string filename = "matrix_" + std::to_string(n) + ".bin";
        binaryFiles.push_back(filename);

        // write performance
        double wr_t = measure(n, true);
        double wr_bw = count_bw(wr_t, matrix_bytes);
        // read performance
        double rd_t = measure(n, false);
        double rd_bw = count_bw(rd_t, matrix_bytes);

        std::cout << std::setw(10) << matrix_bytes / (1024 * 1024) << "MB"
        << std::setw(10) << n << "x" << n 
        << std::setw(10) << std::fixed << std::setprecision(4) << wr_t << "s"
        << std::setw(10) << std::fixed << std::setprecision(2) << wr_bw << "GB/s"
        << std::setw(10) << std::fixed << std::setprecision(4) << rd_t << "s"
        << std::setw(10) << std::fixed << std::setprecision(2) << rd_bw << "GB/s"
        << std::endl;

        csv_file << n << "," 
                 << matrix_mb << "," 
                 << wr_t << "," << wr_bw << "," 
                 << rd_t << "," << rd_bw << "\n";
    }

    csv_file.close();
    clean_files(binaryFiles);
    std::cout << "Results saved. Process completed\n";

    return 0;
}