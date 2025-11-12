#include <iostream>
#include <cmath>
#include <thread>
#include <vector>
#include <mutex>
#include <chrono>
#include <fstream>

double f(double x){
    double df = (1.0 / x) - (x / 4.0);
    double fx = sqrt(1.0 + df*df);
    return fx;
}

void compute_partial_sum(double &sum, std::mutex &sum_mutex, double a, double h, int n, int num_threads, int i){
    double partial_sum = 0;
    int start = (n / num_threads) * i;
    int end = (i == num_threads - 1) ? n : (n / num_threads) * (i + 1);
    for (int j = start; j < end; j++){
        double x = a + j * h;
        partial_sum += f(x);
    }
    partial_sum *= h;

    // lock the mutex and update the sum
    sum_mutex.lock();
    sum += partial_sum;
    sum_mutex.unlock();
}

double parallel_length(double a, double b, int n, int num_threads){
    double h = (b - a) / n;
    double sum = 0.0;
    std::vector<std::thread> threads(num_threads);

    std::mutex sum_mutex;

    for (int i = 0; i < num_threads; ++i){
        threads[i] = std::thread(compute_partial_sum, std::ref(sum), std::ref(sum_mutex), a, h, n, num_threads, i);
    }

    for (int i = 0; i < num_threads; ++i){
        threads[i].join();
    }

    return sum;
}

void test_threads(){
    const double a = 1.0;
    const double b = 6.0;
    const int n = 100000000;
    const std::vector<int> threads_count = {1, 2, 4, 8, 16};
    std::ofstream ofs("threads.csv");
    ofs << "threads,time\n";

    for (int num_threads : threads_count){
        auto start = std::chrono::high_resolution_clock::now();
        parallel_length(a, b, n, num_threads);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        long double elapsed_time = duration.count() * 1e-9;

        ofs << num_threads << "," << elapsed_time << "\n";
        std::cout << "Threads: " << num_threads << ", Time: " << elapsed_time << "s\n";
    }
    std::cout << "Data saved as 'threads.csv' file" << std::endl;
}

void test_numpts(){
    const double a = 1.0;
    const double b = 6.0;
    const int num_threads = 4;
    const std::vector<int> pts = {10, 100, 1000, 10000, 100000, 1000000};
    const double exact = log(6) + 35.0 / 8.0;
    std::ofstream ofs("numpts.csv");
    ofs << "points,error,log_err\n";

    for (int n : pts){
        double approx = parallel_length(a, b, n, num_threads);
        double error = fabs(approx - exact);
        double log_error = log10(error);
        
        ofs << n << "," << error << "," << log_error << "\n";
        std::cout << "n = " << n << ", error = " << error << ", log10(error) = " << log_error << "\n";
    }
    std::cout << "Data saved as 'numpts.csv' file" << std::endl;
}

int main(int argc, char *argv[]){

    const double a = 1.0;
    const double b = 6.0;

    auto start  = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    if (argc != 3){
        std::cerr << "Usage: " << argv[0] << " num_partitions num_threads " << std::endl;
        return 1;
    }
    const int num_partitions = std::atoi(argv[1]);
    const int num_threads = std::atoi(argv[2]);

    start  = std::chrono::high_resolution_clock::now();
    double par_length = parallel_length(a, b, num_partitions, num_threads);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    long double elapsed_time = (duration.count() * 1.e-9);

    std::cout << "Numerical length using C++ threads: " << par_length << "\nTime: " << elapsed_time  << " s" << std::endl;

    std::cout << "_________________________________" << std::endl;
    std::cout << "Test for plotting threads vs time" << std::endl;
    test_threads();
    std::cout << "_________________________________" << std::endl;
    std::cout << "Test for plotting errors vs partitions" << std::endl;
    test_numpts();

    return 0;
}