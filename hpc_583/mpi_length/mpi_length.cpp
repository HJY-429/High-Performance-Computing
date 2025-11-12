#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>

double f(double x){
    double df = (1.0 / x) - (x / 4.0);
    double fx = sqrt(1.0 + df*df);
    return fx;
}

double sequential_sum(double a, double b, int n){
    double sum = 0;
    double h = (b - a) / n;

    for (int i = 0; i < n; ++i){
        double x = a + i * h;
        sum += f(x) * h;
    }
    return sum;
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    int ip, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &ip);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    double a = 1.0, b = 6.0;
    int n = 100000000;

    if (argc < 2) {
        if (ip == 0) {
            std::cerr << "Usage: " << argv[0] << " num_partitions " << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    n = atoi(argv[1]);

    double local_a = a + ip * (b - a) / np;
    double local_b = a + (ip + 1) * (b - a) / np;
    int local_n = n / np;
    double global_sum;

    double start_time = MPI_Wtime();
    double local_sum = sequential_sum(local_a, local_b, local_n);
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    if (ip == 0){
        std::cout << "Numerical length using MPI: " << global_sum << std::endl;
        std::cout << "Time: " << elapsed_time << " s" << std::endl;
        
        // strong scaling efficiency
        if (argc > 2 && std::string(argv[2]) == "--strong-scaling") {
            std::ofstream ofs("strong_scaling_mpi.csv", std::ios_base::app);
            ofs << np << "," << elapsed_time << "\n";
        }
        
        // For numerical error tests
        if (argc > 2 && std::string(argv[2]) == "--numerical-error") {
            const double exact = log(6) + 35.0/8.0;
            double error = fabs(global_sum - exact);
            std::ofstream ofs("numerical_error_mpi.csv", std::ios_base::app);
            ofs << n << "," << error << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
