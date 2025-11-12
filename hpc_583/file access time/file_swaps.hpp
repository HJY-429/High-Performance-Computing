#ifndef FILE_SWAPS_HPP
#define FILE_SWAPS_HPP

#include <fstream>
#include <vector>
#include <stdexcept>

void swapRowsInFile(std::fstream &file, int nRows, int nCols, int i, int j){
    if (i == j) return;
    if (i < 0 || i >= nRows || j < 0 || j >= nRows) {
        throw std::invalid_argument("index i or j must be in between 0 and nRows.");
    }

    std::vector<double> rowI(nCols);
    std::vector<double> rowJ(nCols);

    // read rows i and j
    for (int col = 0; col < nCols; ++col){
        file.seekg((col * nRows + i) * sizeof(double));
        file.read(reinterpret_cast<char *>(&rowI[col]), sizeof(double));
        file.seekg((col * nRows + j) * sizeof(double));
        file.read(reinterpret_cast<char *>(&rowJ[col]), sizeof(double));
    }

    // write swapped rows
    for (int col = 0; col < nCols; ++col){
        file.seekp((col * nRows + i) * sizeof(double));
        // i --> j
        file.write(reinterpret_cast<const char *>(&rowJ[col]), sizeof(double));
        file.seekp((col * nRows + j) * sizeof(double));
        // j --> i
        file.write(reinterpret_cast<const char *>(&rowI[col]), sizeof(double));
    }
}

void swapColsInFile(std::fstream &file, int nRows, int nCols, int i, int j){
    if (i == j) return;
    if (i < 0 || i >= nCols || j < 0 || j >= nCols) {
        throw std::invalid_argument("index i or j must be in between 0 and nCols.");
    }

    std::vector<double> colI(nRows);
    std::vector<double> colJ(nRows);

    // read rows i and j
    file.seekg((i * nRows) * sizeof(double));
    file.read(reinterpret_cast<char *>(colI.data()), nRows * sizeof(double));
    file.seekg((j * nRows) * sizeof(double));
    file.read(reinterpret_cast<char *>(colJ.data()), nRows * sizeof(double));

    // write swapped rows
    file.seekp((i * nRows) * sizeof(double));
    // i --> j
    file.write(reinterpret_cast<const char *>(colJ.data()), nRows * sizeof(double));
    file.seekp((j * nRows) * sizeof(double));
    // j --> i
    file.write(reinterpret_cast<const char *>(colI.data()), nRows * sizeof(double));
}

#endif  
// file_swaps.hpp