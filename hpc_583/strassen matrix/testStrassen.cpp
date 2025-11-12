// hjy429@uw.edu
// AMATH 483-583
// testStrassen.cpp : test performance and correctness of Strassen multiplication

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <cmath>

using namespace std;
using namespace std::chrono;

template <typename T>
vector<vector<T>> addMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B){
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            C[i][j] = A[i][j] + B[i][j];
            }
    }
    return C;
}

template <typename T>
vector<vector<T>> subtractMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B){
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            C[i][j] = A[i][j] - B[i][j];
            }
        }
    return C;
}

template <typename T>
vector<vector<T>> standardMultiply(const vector<vector<T>> &A, const vector<vector<T>> &B) {
    int n = A.size();
    vector<vector<T>> C(n, vector<T>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

size_t nextPowerOfTwo(size_t n) {
    size_t power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

template <typename T>
vector<vector<T>> padMatrix(const vector<vector<T>> &matrix, int newSize) {
    int oldSize = matrix.size();
    vector<vector<T>> padded(newSize, vector<T>(newSize, 0));
    for (int i = 0; i < oldSize; ++i)
        for (int j = 0; j < oldSize; ++j)
            padded[i][j] = matrix[i][j];
    return padded;
}

template <typename T>
vector<vector<T>> unpadMatrix(const vector<vector<T>> &matrix, int oriSize) {
    vector<vector<T>> result(oriSize, vector<T>(oriSize));
    for (int i = 0; i < oriSize; ++i)
        for (int j = 0; j < oriSize; ++j)
            result[i][j] = matrix[i][j];
    return result;
}

template <typename T>
vector<vector<T>> strassenRecursive(const vector<vector<T>> &A, const vector<vector<T>> &B) {
    int n = A.size();
    if (n <= 2) {  // Base case
        return standardMultiply(A, B);
    }

    int k = n / 2;
    
    vector<vector<T>> A11(k, vector<T>(k)), A12(k, vector<T>(k)), A21(k, vector<T>(k)), A22(k, vector<T>(k));
    vector<vector<T>> B11(k, vector<T>(k)), B12(k, vector<T>(k)), B21(k, vector<T>(k)), B22(k, vector<T>(k));

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + k];
            B21[i][j] = B[i + k][j];
            B22[i][j] = B[i + k][j + k];
        }
    }

    auto P1 = strassenRecursive(addMatrix(A11, A22), addMatrix(B11, B22));
    auto P2 = strassenRecursive(addMatrix(A21, A22), B11);
    auto P3 = strassenRecursive(A11, subtractMatrix(B12, B22));
    auto P4 = strassenRecursive(A22, subtractMatrix(B21, B11));
    auto P5 = strassenRecursive(addMatrix(A11, A12), B22);
    auto P6 = strassenRecursive(subtractMatrix(A21, A11), addMatrix(B11, B12));
    auto P7 = strassenRecursive(subtractMatrix(A12, A22), addMatrix(B21, B22));

    auto C11 = addMatrix(subtractMatrix(addMatrix(P1, P4), P5), P7);
    auto C12 = addMatrix(P3, P5);
    auto C21 = addMatrix(P2, P4);
    auto C22 = addMatrix(subtractMatrix(addMatrix(P1, P3), P2), P6);

    vector<vector<T>> C(n, vector<T>(n));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + k] = C12[i][j];
            C[i + k][j] = C21[i][j];
            C[i + k][j + k] = C22[i][j];
        }
    }
    return C;
}

template <typename T>
vector<vector<T>> strassenMultiply(const vector<vector<T>> &A, const vector<vector<T>> &B) {
    int n = A.size();
    int m = nextPowerOfTwo(n);

    auto A_pad = padMatrix(A, m);
    auto B_pad = padMatrix(B, m);
    auto C_pad = strassenRecursive(A_pad, B_pad);

    return unpadMatrix(C_pad, n);
}

template <typename T>
void printMatrix(const vector<vector<T>> &matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

template <typename T>
bool matricesAreEqual(const vector<vector<T>> &A, const vector<vector<T>> &B, double eps = 1e-6)
{
    int n = A.size();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (abs(A[i][j] - B[i][j]) > eps)
                return false;
    return true;
}

template <typename T>
vector<vector<T>> generateRandomMatrix(int n, double start=0.0, double end=1.0) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(start, end);
    
    vector<vector<double>> A(n, vector<double>(n));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = dist(gen);
        }
    }
    return A;
}

void testStrassen(){
    // // dimension = 2, type = <int>
    // vector<vector<int>> A = {{1, 2}, {3, 4}};
    // vector<vector<int>> B = {{5, 6}, {7, 8}};
    // // dimension = 6, type = <double>
    // vector<vector<double>> A = {{1, 2.54, 3, 5, 7, 3}, 
    //                                 {5, 6, 7, 5, 4, 2}, 
    //                                 {9, 10, 11, 3, 1, 2}, 
    //                                 {13, 14, 15, 2, 3, 2}, 
    //                                 {2, 321, 3, 3, 1, 2},  
    //                                 {2, 321, 3, 3, 1, 2}};
    // vector<vector<double>> B = {{16, 15, 14, 2, 1.1, 2}, 
    //                                 {12, 11.111, 102, 2, 1, 1}, 
    //                                 {13, 14, 15, 2, 3, 3}, 
    //                                 {2, 321, 3, 3, 1, 2}, 
    //                                 {2, 4, 1, 2.1, 1, 2},  
    //                                 {2, 321, 3, 3, 1, 2}};

    // test different dimensions
    for (int n = 2; n <= 64; n += 2) {
        auto A = generateRandomMatrix<double>(n, 3, 5.5);
        auto B = generateRandomMatrix<double>(n, 2, 3.21);

        auto C1 = strassenMultiply(A, B);
        auto C2 = standardMultiply(A, B);

        if (n <= 9){
            cout << "______________________________________" << endl;
            cout << "Matrix A:" << endl;
            printMatrix(A);
            cout << "\nMatrix B:" << endl;
            printMatrix(B);
            cout << "\nStandard multiplication result:" << endl;
            printMatrix(C1);
            cout << "\nStrassen multiplication result:" << endl;
            printMatrix(C2);
        }

        if (matricesAreEqual(C1, C2)) {
            cout << "You are GOOD for n = " << n << endl;
        } else {
            cout << "FAILED when n = " << n << endl;
        }
    }
}

void performanceStrassen(){
    const int ntrial = 3;
    const int min_size = 2;
    const int max_size = 64;

    ofstream ofs("performanceStrassen t1.csv");
    ofs << "n,Avg_Time(s),FLOPs\n";
    
    for (int n = min_size; n <= max_size; n += 2) {
        double elapsed_time = 0;
        double avg_time;
        // initialize matrices
        auto A = generateRandomMatrix<double>(n);
        auto B = generateRandomMatrix<double>(n);

        for (int t = 0; t < ntrial; t++) {
            auto start = high_resolution_clock::now();
            
            auto C = strassenMultiply(A, B);
            
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        avg_time = elapsed_time / ntrial;
        // double flops = pow(n, log2(7)) / avg_time;
        double flops = (1.75 * n * n * n + 4.5 * n * n) / avg_time;
        std::cout << "n = " << n 
                  << ", Avg Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << std::endl;
        ofs << n << "," << avg_time << "," << flops << "\n";
    }

    cout << "\nPerformance data saved" << endl;
    ofs.close();
}

int main() {
    testStrassen();
    // performanceStrassen();
    
    return 0;
}

// int
template vector<vector<int>> addMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> subtractMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> standardMultiply<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> strassenMultiply<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template void printMatrix<int>(const vector<vector<int>> &matrix);

// double
template vector<vector<double>> addMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> subtractMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> standardMultiply<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> strassenMultiply<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template void printMatrix<double>(const vector<vector<double>> &matrix);
