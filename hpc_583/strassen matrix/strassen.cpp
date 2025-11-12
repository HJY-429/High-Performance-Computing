// hjy429@uw.edu
// AMATH 483-583
// strassen.cpp : starter code for Strassen implementation

#include <iostream>
#include <vector>

using namespace std;

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
