// hjy429@uw.edu

#ifndef MEM_SWAPS_HPP
#define MEM_SWAPS_HPP

#include <vector>


void swapRows(std::vector<double> &matrix, int nRows, int nCols, int i, int j){
    for (int col = 0; col < nCols; ++col){
        double swap = matrix[col * nRows + i];
        matrix[col * nRows + i] = matrix[col * nRows + j];
        matrix[col * nRows + j] = swap; 
    }
}
void swapCols(std::vector<double> &matrix, int nRows, int nCols, int i, int j){
    for (int row = 0; row < nRows; ++row){
        double swap = matrix[i * nRows + row];
        matrix[i * nRows + row] = matrix[j * nRows + row];
        matrix[j * nRows + row] = swap; 
    }
}

#endif
// MEM_SWAPS_HPP