#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <cstring>
#include <png.h>

extern void grayscaleSequential(png_bytep* image, int width, int height, int channels);

//threaded solution - dummy code
void grayscaleThreaded(png_bytep* image, int width, int height, int channels, int numThreads) {
    std::vector<std::thread> threads;
    int rowsPerThread = height / numThreads;
    int rowsRemain = height % numThreads;
    int startRow = 0;

    for (int i = 0; i < numThreads; ++i){
        int endRow = startRow + rowsPerThread + (i < rowsRemain ? 1 : 0);
        threads.emplace_back([=](){
            png_bytep* subImage = new png_bytep[endRow - startRow];
            for (int y = startRow, idx = 0; y < endRow; y++, idx++){
                subImage[idx] = image[y];
            }
            grayscaleSequential(subImage, width, endRow - startRow, channels);
            delete[] subImage;
        });
        startRow = endRow;
    }

    for (auto& thread : threads){
        thread.join();
    }

}
