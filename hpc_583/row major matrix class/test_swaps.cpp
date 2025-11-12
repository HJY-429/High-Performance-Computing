#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <iomanip>
#include "file_swaps.hpp"

void initMatrix(const std::string& filename, int nRows, int nCols) {
    std::vector<double> matrix(nRows * nCols);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (auto& val : matrix) {
        val = dis(gen);
    }
    
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(matrix.data()), matrix.size() * sizeof(double));
    file.close();
}

void printMatrix(const std::string& filename, int nRows, int nCols) {
    std::ifstream file(filename, std::ios::binary);
    std::vector<double> matrix(nRows * nCols);
    file.read(reinterpret_cast<char*>(matrix.data()), matrix.size() * sizeof(double));
    file.close();

    for (int row = 0; row < nRows; ++row) {
        for (int col = 0; col < nCols; ++col) {
            std::cout << std::setw(8) << matrix[col * nRows + row] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int main() {
    const int ntrial = 3;
    // const std::vector<int> sizes = {3};      // test swap
    const std::vector<int> sizes = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};

    std::ofstream csvFile("swap_matrix.csv");
    csvFile << "Size,AvgRowSwapTime,AvgColSwapTime\n";
    std::cout << "Size,AvgRowSwapTime,AvgColSwapTime\n";
    
    for (int n : sizes) {
        const std::string filename = "matrix_" + std::to_string(n) + ".bin";
        double totalRowTime = 0.0;
        double totalColTime = 0.0;
        
        for (int trial = 0; trial < ntrial; ++trial) {
            initMatrix(filename, n, n);

            if (n <= 5) {
                std::cout << "  Original Matrix (" << n << "x" << n << ", Trial " << trial + 1 << ")  \n";
                printMatrix(filename, n, n);
            }
            
            // Test row swap
            {
            std::fstream file(filename, std::ios::in | std::ios::out | std::ios::binary);
            auto startRow = std::chrono::high_resolution_clock::now();
            swapRowsInFile(file, n, n, 0, n-1); // 1 <--> n
            auto endRow = std::chrono::high_resolution_clock::now();
            totalRowTime += std::chrono::duration<double>(endRow - startRow).count();
            file.close();

            if (n <= 5) {
                std::cout << "  After Row Swap  \n";
                printMatrix(filename, n, n);
            }
            }
            
            // Test column swap
            {
            std::fstream file(filename, std::ios::in | std::ios::out | std::ios::binary);
            auto start = std::chrono::high_resolution_clock::now();
            swapColsInFile(file, n, n, 0, n-1); // 1 <--> n
            auto end = std::chrono::high_resolution_clock::now();
            totalColTime += std::chrono::duration<double>(end - start).count();
            file.close();

            if (n <= 5) {
                std::cout << "  After Column Swap  \n";
                printMatrix(filename, n, n);
            }
            }
            
            std::remove(filename.c_str());
        }
        
        double avgRowTime = totalRowTime / ntrial;
        double avgColTime = totalColTime / ntrial;
        csvFile << n << "," << avgRowTime << "," << avgColTime << "\n";
        std::cout << n << "," << avgRowTime << "," << avgColTime << "\n";
    }
    csvFile.close();
    
    return 0;
}