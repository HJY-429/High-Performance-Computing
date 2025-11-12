Compilation Instructions for librefBLAS.so:

(Make sure you are in the same path of this file: README.txt)

0. Clean previous builds
   rm *.o *.so

1. First compile each source file to position-independent object code:
   g++ -std=c++17 -c -fPIC ref_daxpy.cpp -o ref_daxpy.o
   g++ -std=c++17 -c -fPIC ref_dgemv.cpp -o ref_dgemv.o
   g++ -std=c++17 -c -fPIC ref_dgemm.cpp -o ref_dgemm.o
   g++ -std=c++17 -c -fPIC ref_axpyt.cpp -o ref_axpyt.o
   g++ -std=c++17 -c -fPIC ref_gemvt.cpp -o ref_gemvt.o
   g++ -std=c++17 -c -fPIC ref_gemmt.cpp -o ref_gemmt.o

2. Combine object files into shared library:
   g++ -shared -o librefBLAS.so ref_daxpy.o ref_dgemv.o ref_dgemm.o ref_axpyt.o ref_gemvt.o ref_gemmt.o

3. Compile the test program: (assuming test program is test.cpp)
   g++ -std=c++17 -o xtest test.cpp -L. -lrefBLAS

4. To run the test program (may need to set LD_LIBRARY_PATH):
   export LD_LIBRARY_PATH=. :$LD_LIBRARY_PATH
   env | grep LD_LIBRARY_PATH
   ./xtest






