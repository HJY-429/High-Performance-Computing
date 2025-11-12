#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

// JKI loop order
template <typename T>
void mm_jki(T a, const std::vector<T>& A, const std::vector<T>& B, T b,
        std::vector<T>& C, int m, int p, int n){
    if (b != 1.0) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                C[i*n + j] *= b;
            }
        }
    }
    
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < p; ++k) {
            T temp = a * B[k*n + j];
            for (int i = 0; i < m; ++i) {
                C[i*n + j] += temp * A[i*p + k];
            }
        }
    }
}
// KIJ loop order
template <typename T>
void mm_kij(T a, const std::vector<T>& A, const std::vector<T>& B, T b,
            std::vector<T>& C, int m, int p, int n) {
    if (b != 1.0) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                C[i*n + j] *= b;
            }
        }
    }
    
    for (int k = 0; k < p; ++k) {
        for (int i = 0; i < m; ++i) {
            T temp = a * A[i*p + k];
            for (int j = 0; j < n; ++j) {
                C[i*n + j] += temp * B[k*n + j];
            }
        }
    }
}

template <typename T>
void measure_jki_perform(const std::string& filename, int ntrials, int max_n) {
    std::ofstream outfile(filename);
    outfile << "n,avg_time,flops\n";
    
    for (int n = 2; n <= max_n; n++) {
        int m = n, p = n;
        std::vector<T> A(m*p, 2.0);
        std::vector<T> B(p*n, 3.0);
        std::vector<T> C(m*n, 1.0);
        T alpha = 2.0, beta = 3.0;

        long double elapsed_time = 0.L;
        
        for (int t = 0; t < ntrials; t++) {
            std::fill(C.begin(), C.end(), 0.0);
            
            auto start = std::chrono::high_resolution_clock::now();
            mm_jki(alpha, A, B, beta, C, m, p, n);
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        
        long double avg_time = elapsed_time / static_cast<long double>(ntrials);
        double flops = (2.0 * n * n * n + 2.0 * n * n) / avg_time;
        
        std::cout << "JKI n = " << n 
                  << ", Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << "\n";
                  
        outfile << n << "," << avg_time << "," << flops << "\n";
    }
    outfile.close();
}

template <typename T>
void measure_kij_perform(const std::string& filename, int ntrials, int max_n) {
    std::ofstream outfile(filename);
    outfile << "n,avg_time,flops\n";
    
    for (int n = 2; n <= max_n; n++) {
        int m = n, p = n;
        std::vector<T> A(m*p, 2.0);
        std::vector<T> B(p*n, 3.0);
        std::vector<T> C(m*n, 1.0);
        T alpha = 2.0, beta = 3.0;

        long double elapsed_time = 0.L;
        
        for (int t = 0; t < ntrials; t++) {
            std::fill(C.begin(), C.end(), 0.0);
            
            auto start = std::chrono::high_resolution_clock::now();
            mm_kij(alpha, A, B, beta, C, m, p, n);
            auto stop = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9);
        }
        
        long double avg_time = elapsed_time / static_cast<long double>(ntrials);
        double flops = (2.0 * n * n * n + 2.0 * n * n) / avg_time;
        
        std::cout << "KIJ n = " << n 
                  << ", Time = " << avg_time << " s"
                  << ", FLOPs = " << flops << "\n";
                  
        outfile << n << "," << avg_time << "," << flops << "\n";
    }
    outfile.close();
}

template <typename T>
void test_mm(int m, int p, int n) {
    std::vector<T> A(m*p, 2.0);
    std::vector<T> B(p*n, 3.0);
    std::vector<T> C_orig(m*n, 1.0);
    
    T alpha = 2.0, beta = 3.0;

    std::cout << "\nValidation for " << (sizeof(T) == 4 ? "float" : "double") 
              << " (" << m << "x" << p << ") * (" << p << "x" << n << ")\n";

    // Reference
    std::vector<T> C_ref = C_orig;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            T sum = 0;
            for (int k = 0; k < p; ++k){
                sum += A[i*p + k] * B[k*n + j];
            }
            C_ref[i*n + j] = alpha * sum + beta * C_ref[i*n + j];
        }
    }

    // Test JKI
    std::vector<T> C_jki = C_orig;
    mm_jki(alpha, A, B, beta, C_jki, m, p, n);
    
    // Test KIJ
    std::vector<T> C_kij = C_orig;
    mm_kij(alpha, A, B, beta, C_kij, m, p, n);

    std::cout << "Reference\tJKI\t\tKIJ\n";
    for (int i = 0; i < m; ++i) {  
        for (int j = 0; j < n; ++j) {
            std::cout << C_ref[i*n + j] << "\t"
                      << C_jki[i*n + j] << "\t"
                      << C_kij[i*n + j] << "\n";
        }
    }
}

int main() {
    const int ntrials = 3;
    const int max_n = 512;

    // test_mm<float>(3, 3, 3);
    // test_mm<double>(3, 3, 3);
    // measure_kij_perform<float>("results_kij_floatO0.csv", ntrials, max_n);
    // measure_jki_perform<float>("results_jki_floatO0.csv", ntrials, max_n);
    // measure_kij_perform<double>("results_kij_doubleO0.csv", ntrials, max_n);
    measure_jki_perform<double>("results_jki_doubleO0.csv", ntrials, max_n);
    
    return 0;
}



// template <typename T>
// void measure_perform(const std::string& filename, int ntrials, int max_n) {
//     std::ofstream outfile(filename);
//     outfile << "n,avg_time_jki,flops_jki,avg_time_kij,flops_kij\n";

//     for (int n = 2; n <= max_n; n++) {
//         int m = n, p = n;
//         std::vector<T> A(m*p, 2.0);
//         std::vector<T> B(p*n, 3.0);
//         std::vector<T> C(m*n, 1.0);
//         T alpha = 2.0, beta = 3.0;

//         long double avg_time_jki;
//         long double avg_time_kij;
//         long double elapsed_time_jki = 0.L;
//         long double elapsed_time_kij = 0.L;
        
//         for (int t = 0; t < ntrials; t++) {

//             std::fill(C.begin(), C.end(), 0.0);
//             // Measure jki
//             auto start = std::chrono::high_resolution_clock::now();
//             mm_jki(alpha, A, B, beta, C, m, p, n);
//             auto stop = std::chrono::high_resolution_clock::now();
//             auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
//             elapsed_time_jki += (duration.count() * 1.e-9);
            
//             // Reset C
//             std::fill(C.begin(), C.end(), 0.0);
//             // Measure kij
//             start = std::chrono::high_resolution_clock::now();
//             mm_kij(alpha, A, B, beta, C, m, p, n);
//             stop = std::chrono::high_resolution_clock::now();
//             duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
//             elapsed_time_kij += (duration.count() * 1.e-9);
//         }
        
//         // FLOPs
//         avg_time_jki = elapsed_time_jki / static_cast<long double>(ntrials);
//         avg_time_kij = elapsed_time_kij / static_cast<long double>(ntrials);
//         double flops_jki = (2.0 * n * n * n + 2.0 * n * n) / avg_time_jki;
//         double flops_kij = (2.0 * n * n * n + 2.0 * n * n) / avg_time_kij;
        
//         // Output results
//         std::cout << "n = " << n
//                   << ", Avg Time JKI = " << avg_time_jki << " s"
//                   << ", FLOPs JKI = " << flops_jki
//                   << ", Avg Time KIJ= " << avg_time_kij << " s"
//                   << ", FLOPs KIJ= " << flops_kij << std::endl;
//         outfile << n << "," 
//         << avg_time_jki << "," << flops_jki << ","
//         << avg_time_kij << "," << flops_kij << "\n";
//     }
//     outfile.close();
// }