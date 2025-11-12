# Repository Subfolders Summary: AMATH 583 - High Performance Computing

## Overview
This repository contains a comprehensive collection of projects and assignments for High Performance Computing (HPC) that covers parallel programming paradigms, performance analysis, and optimization techniques. The projects span multi-threading, distributed computing, GPU programming, and hybrid approaches.

---

## Subfolder Descriptions and Problems Solved

### 1. **Threaded Image Filter** (`threaded image filter/`)
**Goal**: Efficient parallel processing of PNG images using multi-threading for grayscale conversion.

**Approach**:
- Implements parallel PNG image filtering using pthreads
- Distributes image rows among threads for scalable processing
- Provides correctness verification through pixel-by-pixel comparison
- Demonstrates performance scaling from 1 to 4 threads

**Key Files**:
- `cpp-png-gs.cpp`: Main sequential implementation and orchestration
- `grayscaleThreaded.cpp`: Multi-threaded implementation
- `compare_png.cpp`: Correctness validation utility
- Performance data showing 1.8x speedup with 4 threads

### 2. **Shared Object Library** (`shared_object_lib/`) 
**Goal**: Development of reference BLAS (Basic Linear Algebra Subprograms) implementations for performance comparison and educational purposes.

**Approach**:
- Complete BLAS implementation as dynamically loadable library (`librefBLAS.so`)
- Threaded and non-threaded variants of core operations
- Systematic performance benchmarking across problem sizes

**Key Components**:
- BLAS Level 1: `ref_daxpy.cpp` (vector-vector operations)
- BLAS Level 2: `ref_dgemv.cpp` (matrix-vector multiplication)
- BLAS Level 3: `ref_dgemm.cpp` (matrix-matrix multiplication)
- Threaded variants: `ref_axpyt.cpp`, `ref_gemvt.cpp`, `ref_gemmt.cpp`
- Performance measurement infrastructure with CSV output

### 3. **Strong Scaling Study** (`strong scaling for parallel solver/`)
**Goal**: Analysis of parallel solver scalability for distributed linear systems.

**Approach**:
- MPI-based strong scaling analysis using ScaLAPACK
- Tests problem sizes N=4096 and N=8192 with 1-16 processes
- Grid-based process topology optimization
- Performance measurement on HYAK cluster

**Key Implementation**:
- `kr-azeb.c`: Parallel linear solver implementation
- Systematic resource allocation with salloc
- Performance tracking showing near-linear scaling for N=4096

### 4. **Performance Comparison Studies**

#### **Compare OpenBLAS to CUBLAS** (`Compare OpenBLAS to CUBLAS on HYAK/`)
**Goal**: CPU vs GPU performance comparison for linear algebra operations.

**Approach**: Benchmarks OpenBLAS (CPU) vs CUBLAS (GPU) libraries with systematic performance measurement and visualization.

#### **Compare FFTW to CUFFT** (`Compare FFTW to CUFFT on HYAK/`)
**Goal**: Fast Fourier Transform performance analysis comparing CPU and GPU implementations.

**Approach**: Direct comparison of FFTW (CPU) and CUFFT (GPU) libraries across different transform sizes.

#### **CPU-GPU Data Copy Speed** (`CPU-GPU data copy speed on HYAK/`)
**Goal**: Memory bandwidth analysis between CPU and GPU memory spaces.

**Approach**: Measures and visualizes data transfer rates between host and device memory.

### 5. **Parallel Algorithms and Implementations**

#### **Strassen Matrix Multiplication** (`strassen matrix/`)
**Goal**: Implementation of fast matrix multiplication algorithm with better asymptotic complexity.

**Approach**: Strassen's algorithm implementation with performance comparison against standard methods.

#### **Templated GEMM** (`templated gemm/`)
**Goal**: Generic matrix multiplication using C++ templates. Compiler Optimization of Matrix Multiplication Loop Permutations.

**Approach**: Template-based implementation for type-generic matrix operations, applying compiler optimization levels -O0 and -O3.

#### **Row Major Matrix Class** (`row major matrix class/`)
**Goal**: Object-oriented matrix abstraction with optimized memory layout.

**Approach**: Custom matrix class implementation optimized for row-major access patterns.

### 6. **System Performance Analysis**

#### **Memory Access Time** (`memory access time/`)
**Goal**: Perform row and column swap operations on a type double matrix stored in column major index order using a single vector container for the data.

**Approach**: Conduct a performance test for square matrix dimensions 16, 32, 64, ... 4096, measuring the time required to conduct row swaps and column swaps separately.

#### **IO Bandwidth** (`IO bandwidth/`)
**Goal**: Explore input/output system performance.

**Approach**: Measures and visualizes I/O bandwidth capabilities.

#### **File Access Time** (`file access time/`)
**Goal**: File system performance analysis.

**Approach**: Benchmarks file access patterns and throughput. Conduct a performance test for square matrix dimensions 16, 32, 64, 128, ... 8192, measuring the time required to conduct file-based row and column swaps separately.

#### **Communication Bandwidth** (`communication bandwidth/`)
**Goal**: Perform a MPI C++ function using point-to-point communication that implements broadcast.

**Approach**: Measures MPI communication bandwidth and latency.

### 7. **Parallel Programming**

#### **Elevator Scheduling** (`elevator/`)
**Goal**: Concurrent programming with resource allocation and synchronization. Perform a threaded elevator scheduler for a 50-story business building with 6 elevators.

**Approach**: Multi-threaded elevator scheduling algorithm implementation.

#### **MPI Length Calculation** (`mpi_length/`)
**Goal**: Distributed computation using message passing.

**Approach**: MPI-based implementation for collaborative length calculations.

#### **Threads Length** (`threads_length/`)
**Goal**: Multi-threaded computation with shared memory.

**Approach**: Thread-based implementation of length calculation algorithms.

#### **OpenBLAS** (`openblas/`)
**Goal**: Utilization of optimized BLAS libraries. Measure and plot the performance of double precision L1, L2, and L3 BLAS using the OpenBLAS library on Hyak.

**Approach**: OpenBLAS library usage examples and performance testing.

### 8. **LAPACK Linear System Solver** (`LAPACK linear system solver/`)
**Goal**: Numerical solution of linear systems using high-performance libraries.

**Approach**: LAPACK-based implementation with performance analysis and optimization studies.

---

## Common Methodology and Performance Analysis

All subfolders demonstrate consistent approaches to:

### **Performance Measurement**:
- High-resolution timing infrastructure
- Multiple trial averaging for statistical significance
- Both wall time and CPU time measurements
- FLOPS (Floating Point Operations Per Second) calculation

### **Scalability Analysis**:
- Strong scaling studies (fixed problem size, varying resources)
- Weak scaling considerations
- Performance visualization using PNG plots
- Identification of scaling limitations and bottlenecks

### **Correctness Validation**:
- Comparison between sequential and parallel implementations
- Numerical accuracy verification
- Application-specific validation (pixel comparison for images, residual checking for numerical methods)

### **Resource Utilization**:
- CPU utilization analysis
- GPU memory management and optimization
- MPI process topology optimization
- Thread-level parallelization strategies