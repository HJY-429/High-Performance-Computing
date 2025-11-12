#include <iostream>

int main(){
    int a = 200;
    int b = 300;
    int c = 400;
    int d = 500;
    int result = a * b * c * d;

    std::cout << "Q3 Result: " << result << std::endl;


    unsigned int counter = 0;
    for (int i = 0; i < 3; ++i){
        --counter;
    }
    std::cout << "Q4 Result: " << counter << std::endl;

    return 0;
}