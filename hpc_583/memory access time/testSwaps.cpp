#include <iostream>
#include <utility>
#include <chrono>
#include <random>
#include <fstream>
#include "mem_swaps.hpp"

std::pair<int, int> getRandomIndices(int n){
    int i = std::rand() % n;
    int j = std::rand() % (n - 1);
    if (j >= i){
        j++;
    }
    return std::make_pair(i, j);
}

void printMatrix(const std::vector<double>& matrix, int nRows, int nCols) {
    for (int row = 0; row < nRows; ++row) {
        for (int col = 0; col < nCols; ++col) {
            std::cout << std::setw(8) << matrix[col * nRows + row] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void testSwaps(){
    const int nRows = 4;
    const int nCols = 5;

    std::vector<double> matrix(nRows * nCols);
    for (int col = 0; col < nCols; ++col) {
        for (int row = 0; row < nRows; ++row) {
            matrix[col * nRows + row] = row * 10.0 + col;
        }
    }
    // test row swaps
    auto rowI = getRandomIndices(nRows);
    std::cout << "Before row swap:\n";
    printMatrix(matrix, nRows, nCols);
    swapRows(matrix, nRows, nCols, rowI.first, rowI.second);
    std::cout << "After row swap (" << rowI.first << "<->" << rowI.second << ") :\n";
    printMatrix(matrix, nRows, nCols);

    // test col swaps
    auto colI = getRandomIndices(nCols);
    std::cout << "Before column swap:\n";
    printMatrix(matrix, nRows, nCols);
    swapCols(matrix, nRows, nCols, colI.first, colI.second);
    std::cout << "After column swap (" << colI.first << "<->" << colI.second << ") :\n";
    printMatrix(matrix, nRows, nCols);
}

void performanceSwaps(){
    const int ntrial = 3;
    std::vector<int> sizes = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};

    std::ofstream ofs("swap_time t1.csv");
    ofs << "MatrixSize,RowSwapTime,ColumnSwapTime\n";
    for (int n : sizes) {
        std::cout << "Start performance of size: " << n << std::endl;
        std::vector<double> matrix(n * n);

        // initialize
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.0, 1.0);
        for (auto& val : matrix) {
            val = dist(gen);
        }

        // row swaps time
        auto rowI = getRandomIndices(n);
        auto startRow = std::chrono::high_resolution_clock::now();
        for (int t = 0; t < ntrial; ++t) {
            swapRows(matrix, n, n, rowI.first, rowI.second);
        }
        auto endRow = std::chrono::high_resolution_clock::now();
        double rowTime = std::chrono::duration<double>(endRow - startRow).count() / ntrial;
        // col swaps time
        auto colI = getRandomIndices(n);
        auto startCol = std::chrono::high_resolution_clock::now();
        for (int t = 0; t < ntrial; ++t) {
            swapCols(matrix, n, n, colI.first, colI.second);
        }
        auto endCol = std::chrono::high_resolution_clock::now();
        double colTime = std::chrono::duration<double>(endCol - startCol).count() / ntrial;

        // write csv files
        ofs << n << "," << rowTime << "," << colTime << "\n";
        std::cout << "Completed size: " << n << std::endl;
    }
}

int main(){
    testSwaps();
    // performanceSwaps();
    return 0;
}