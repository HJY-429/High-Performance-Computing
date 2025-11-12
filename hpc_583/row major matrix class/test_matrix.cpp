#include <iostream>
#include "matrix_class.hpp"

int main() {
    std::cout << "Testing Matrix Transpose:\n";
    Matrix<int> A(2, 3);
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 4; A(1, 1) = 5; A(1, 2) = 6;
    
    Matrix<int> AT = A.transpose();
    std::cout << "Matrix A:\n";
    for (int i = 0; i < A.numRows(); ++i) {
        for (int j = 0; j < A.numCols(); ++j) {
            std::cout << A(i, j) << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Transposed Matrix A^T:\n";
    for (int i = 0; i < AT.numRows(); ++i) {
        for (int j = 0; j < AT.numCols(); ++j) {
            std::cout << AT(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";


    std::cout << "Testing Infinity Norm:\n";
    Matrix<float> B(2, 2);
    B(0, 0) = 1.0f; B(0, 1) = -2.0f;
    B(1, 0) = -3.30f; B(1, 1) = 4.10f;
    
    double norm = B.infinityNorm();
    std::cout << "Matrix:\n";
    for (int i = 0; i < B.numRows(); ++i) {
        for (int j = 0; j < B.numCols(); ++j) {
            std::cout << B(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Infinity Norm: " << norm << "\n\n";


    std::cout << "Testing Matrix Multiplication:\n";
    Matrix<int> C(2, 3);
    C(0, 0) = 1; C(0, 1) = 2; C(0, 2) = 3;
    C(1, 0) = 4; C(1, 1) = 5; C(1, 2) = 6;
    
    Matrix<int> D(3, 2);
    D(0, 0) = 7; D(0, 1) = 8;
    D(1, 0) = 9; D(1, 1) = 10;
    D(2, 0) = 11; D(2, 1) = 12;
    
    Matrix<int> CD = C * D;
    std::cout << "Matrix C:\n";
    for (int i = 0; i < C.numRows(); ++i) {
        for (int j = 0; j < C.numCols(); ++j) {
            std::cout << C(i, j) << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Matrix D:\n";
    for (int i = 0; i < D.numRows(); ++i) {
        for (int j = 0; j < D.numCols(); ++j) {
            std::cout << D(i, j) << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Matrix C*D:\n";
    for (int i = 0; i < CD.numRows(); ++i) {
        for (int j = 0; j < CD.numCols(); ++j) {
            std::cout << CD(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    
    std::cout << "Testing Matrix Addition:\n";
    Matrix<double> E(2, 2);
    E(0, 0) = 1.2; E(0, 1) = 2;
    E(1, 0) = 3; E(1, 1) = 4;
    
    Matrix<double> F(2, 2);
    F(0, 0) = 5.23; F(0, 1) = 6;
    F(1, 0) = 7; F(1, 1) = 8.12;
    
    Matrix<double> EF = E + F;
    std::cout << "Matrix E:\n";
    for (int i = 0; i < E.numRows(); ++i) {
        for (int j = 0; j < E.numCols(); ++j) {
            std::cout << E(i, j) << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Matrix F:\n";
    for (int i = 0; i < F.numRows(); ++i) {
        for (int j = 0; j < F.numCols(); ++j) {
            std::cout << F(i, j) << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "Matrix E+F:\n";
    for (int i = 0; i < EF.numRows(); ++i) {
        for (int j = 0; j < EF.numCols(); ++j) {
            std::cout << EF(i, j) << " ";
        }
        std::cout << "\n";
    }

    return 0;
}