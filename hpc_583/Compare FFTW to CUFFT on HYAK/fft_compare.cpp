#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <fftw3.h>
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>

const double PI = 3.14159265359;

void init_k_vec(std::vector<double>&kx, std::vector<double>&ky, std::vector<double>&kz,
    int nx, int ny, int nz, double lx, double ly, double lz){
        std::vector<double> kkx(nx), kky(ny), kkz(nz);
        int i = 0;
        
        // init kx components
        kkx[0] = 0.0;
        for (int i = 1; i <= nx / 2 - 1; i++){
            kkx[i] = 2.0 * PI / lx * static_cast<double>(i);
        }
        int j = -i;
        for (int i = nx / 2; i < nx; i++){
            kkx[i] = 2.0 * PI / lx * static_cast<double>(j);
            j++;
        }

        // init ky components
        kky[0] = 0.0;
        for (int i = 1; i <= ny / 2 - 1; i++){
            kky[i] = 2.0 * PI / ly * static_cast<double>(i);
        }
        j = -i;
        for (int i = ny / 2; i < ny; i++){
            kky[i] = 2.0 * PI / ly * static_cast<double>(j);
            j++;
        }

        // init kz components
        kkz[0] = 0.0;
        for (int i = 1; i <= nz / 2 - 1; i++){
            kkz[i] = 2.0 * PI / lz * static_cast<double>(i);
        }
        j = -i;
        for (int i = nz / 2; i < nz; i++){
            kkz[i] = 2.0 * PI / lz * static_cast<double>(j);
            j++;
        }

        // fill 3D k-vectors
        for (int ix = 0; ix < nx; ix++){
            for (int iy = 0; iy < ny; iy++){
                for (int iz = 0; iz < nz; iz++){
                    int idx = iz + nz * (iy + ny * ix);
                    kx[idx] = kkx[ix];
                    ky[idx] = kky[iy];
                    kz[idx] = kkz[iz];
                }
            }
        }

}
void init_wave(std::vector<std::complex<double>>& wave, const std::vector<double>&x, const std::vector<double>&y,
    const std::vector<double>&z, int nx, int ny, int nz, double kx, double ky, double kz) {
        int nxyz = nx * ny * nz;
        for (int i = 0; i < nxyz; i++){
            wave[i] = std::exp(std::complex<double>(0.0, 1.0) * (x[i] * kx + y[i] * ky + z[i] * kz));
        }

        double norm = 0.0;
        for (int i = 0; i < nxyz; i++){
            norm += std::norm(wave[i]);
        }
        double sqrt_norm = std::sqrt(norm);
        for (int i = 0; i < nxyz; i++){
            wave[i] /= sqrt_norm;
        }
}


double FFTW_gradient(int nx, int ny, int nz, int ntrials){
    int nxyz = nx * ny * nz;
    double lx = 1.0, ly = 1.0, lz = 1.0;
    std::vector<double> kx(nxyz), ky(nxyz), kz(nxyz);
    init_k_vec(kx, ky, kz, nx, ny, nz, lx, ly, lz);
    std::vector<double> x(nxyz), y(nxyz), z(nxyz);
    double center_x = static_cast<double>(nx) / 2.0;
    double center_y = static_cast<double>(ny) / 2.0;
    double center_z = static_cast<double>(nz) / 2.0;
    double dx = lx / nx, dy = ly / ny, dz = lz / nz;

    for (int ix = 0; ix < nx; ix++){
        for (int iy = 0; iy < ny; iy++){
            for (int iz = 0; iz < nz; iz++){
                int idx = iz + nz * (iy + ny * ix);
                x[idx] = dx * static_cast<double>(ix) - center_x * dx;
                y[idx] = dy * static_cast<double>(iy) - center_y * dy;
                z[idx] = dz * static_cast<double>(iz) - center_z * dz;
            }
        }
    }

    std::vector<std::complex<double>> wave(nxyz);
    init_wave(wave, x, y, z, nx, ny, nz, 2.0 * PI, 2.0 * PI, 2.0 * PI);

    std::vector<std::complex<double>> fft_3(nxyz), d_tmp(nxyz);
    std::vector<std::complex<double>> d_dx(nxyz), d_dy(nxyz), d_dz(nxyz);
    fftw_plan planF = fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(&wave[0]), 
    reinterpret_cast<fftw_complex*>(&fft_3[0]), FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan planB = fftw_plan_dft_3d(nx, ny, nz, reinterpret_cast<fftw_complex*>(&fft_3[0]), 
    reinterpret_cast<fftw_complex*>(&d_tmp[0]), FFTW_BACKWARD, FFTW_MEASURE);
    
    double elapsed_time = 0.0;
    double avg_time;

    for (int t = 0; t < ntrials; t++) {
        auto start = std::chrono::high_resolution_clock::now();
        
        fftw_execute(planF);

        for (int j = 0; j < nxyz; j++){
            fft_3[j] = wave[j] / static_cast<double>(nxyz);
            fft_3[j] = fft_3[j] * std::complex<double>(0.0, kx[j]);
        }
        fftw_execute(planB);
        for (int j = 0; j < nxyz; j++){
            d_dx[j] = d_tmp[j] / static_cast<double>(nxyz);
        }

        for (int j = 0; j < nxyz; j++){
            fft_3[j] = fft_3[j] * std::complex<double>(0.0, ky[j]);
        }
        fftw_execute(planB);
        for (int j = 0; j < nxyz; j++){
            d_dy[j] = d_tmp[j] / static_cast<double>(nxyz);
        }

        for (int j = 0; j < nxyz; j++){
            fft_3[j] = fft_3[j] * std::complex<double>(0.0, kz[j]);
        }
        fftw_execute(planB);
        for (int j = 0; j < nxyz; j++){
            d_dz[j] = d_tmp[j] / static_cast<double>(nxyz);
        }
        
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        elapsed_time += (duration.count() * 1.e-9);
    }

    // delete
    fftw_destroy_plan(planF);
    fftw_destroy_plan(planB);
    
    avg_time = elapsed_time / ntrials;
    double flops = (24 * nxyz * log2(nxyz) + 9 * nxyz) / avg_time;
    
    return flops;
}

__global__ void scale_kernel(cufftDoubleComplex* data, int n, double* k, int offset, double scale){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n){
        double k_val = k[offset * n + idx];
        data[idx].x *= -k_val * scale;
        data[idx].y *= k_val * scale;
    }
}
double CUFFT_gradient(int nx, int ny, int nz, int ntrials){
    int nxyz = nx * ny * nz;
    double lx = 1.0, ly = 1.0, lz = 1.0;
    std::vector<double> kx(nxyz), ky(nxyz), kz(nxyz);
    init_k_vec(kx, ky, kz, nx, ny, nz, lx, ly, lz);

    std::vector<double> x(nxyz), y(nxyz), z(nxyz);
    double center_x = static_cast<double>(nx) / 2.0;
    double center_y = static_cast<double>(ny) / 2.0;
    double center_z = static_cast<double>(nz) / 2.0;
    double dx = lx / nx, dy = ly / ny, dz = lz / nz;

    for (int ix = 0; ix < nx; ix++){
        for (int iy = 0; iy < ny; iy++){
            for (int iz = 0; iz < nz; iz++){
                int idx = iz + nz * (iy + ny * ix);
                x[idx] = dx * static_cast<double>(ix) - center_x * dx;
                y[idx] = dy * static_cast<double>(iy) - center_y * dy;
                z[idx] = dz * static_cast<double>(iz) - center_z * dz;
            }
        }
    }

    std::vector<std::complex<double>> wave(nxyz);
    init_wave(wave, x, y, z, nx, ny, nz, 2.0 * PI, 2.0 * PI, 2.0 * PI);
    cufftDoubleComplex *d_wave, *d_fft, *d_dx, *d_dy, *d_dz;
    double *d_kxyz;

    cudaMalloc((void**)&d_wave, sizeof(cufftDoubleComplex) * nxyz);
    cudaMalloc((void**)&d_fft, sizeof(cufftDoubleComplex) * nxyz);
    cudaMalloc((void**)&d_dx, sizeof(cufftDoubleComplex) * nxyz);
    cudaMalloc((void**)&d_dy, sizeof(cufftDoubleComplex) * nxyz);
    cudaMalloc((void**)&d_dz, sizeof(cufftDoubleComplex) * nxyz);
    cudaMalloc((void**)&d_kxyz, sizeof(double) * 3 * nxyz);

    cudaMemcpy(d_wave, wave.data(), sizeof(cufftDoubleComplex) * nxyz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_kxyz, kx.data(), sizeof(double) * nxyz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_kxyz + nxyz, ky.data(), sizeof(double) * nxyz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_kxyz + 2 * nxyz, kz.data(), sizeof(double) * nxyz, cudaMemcpyHostToDevice);

    cufftHandle plan;
    cufftPlan3d(&plan, nx, ny, nz, CUFFT_Z2Z);

    int threadsPerBlock = 256;
    int blocks = (nxyz + threadsPerBlock - 1) / threadsPerBlock;

    double elapsed_time = 0.0;
    double avg_time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    for (int t = 0; t < ntrials; t++) {
        cudaEventRecord(start);
        
        cufftExecZ2Z(plan, d_wave, d_fft, CUFFT_FORWARD);

        scale_kernel<<<blocks, threadsPerBlock>>>(d_fft, nxyz, d_kxyz, 0, 1.0 / nxyz);
        cufftExecZ2Z(plan, d_fft, d_dx, CUFFT_INVERSE);

        scale_kernel<<<blocks, threadsPerBlock>>>(d_fft, nxyz, d_kxyz, 1, 1.0 / nxyz);
        cufftExecZ2Z(plan, d_fft, d_dy, CUFFT_INVERSE);

        scale_kernel<<<blocks, threadsPerBlock>>>(d_fft, nxyz, d_kxyz, 2, 1.0 / nxyz);
        cufftExecZ2Z(plan, d_fft, d_dz, CUFFT_INVERSE);
        
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        elapsed_time += milliseconds;
    }

    // delete
    cufftDestroy(plan);
    cudaFree(d_wave);
    cudaFree(d_fft);
    cudaFree(d_dx);
    cudaFree(d_dy);
    cudaFree(d_dz);
    cudaFree(d_kxyz);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    avg_time = (elapsed_time / ntrials) * 1.e-3;
    double flops = (24 * nxyz * log2(nxyz) + 9 * nxyz) / avg_time;

    return flops;
}

int main(){
    const int ntrials = 3;
    const int min_n = 16;
    const int max_n = 256;

    // Create and open a CSV file
    std::ofstream outfile("performance.csv");
    outfile << "n,FFTW_FLOPs,CUFFT_FLOPs\n";

    for (int n = min_n; n <= max_n; n *= 2){
        std::cout << "Running for n = " << n << std::endl;

        double fftw_flops = FFTW_gradient(n, n, n, ntrials);
        double cufft_flops = CUFFT_gradient(n, n, n, ntrials);

        std::cout << "n: " << n << ", FFTW: " << fftw_flops << ", CUFFT: " << cufft_flops << std::endl;
        outfile << n << "," << fftw_flops << "," << cufft_flops << "\n";
    }

    outfile.close();

    return 0;
}
