#include <iostream>
#include <iomanip>
#include <cmath>

int main(){
    // SP (32 bit)
    float eps_f = 1.0f;
    int j1 = 0;
    while ((1.0f + eps_f / 2.0f) != 1.0f){
        eps_f /= 2.0f;
        j1++;
    }

    // DP (64 bit)
    double eps_d = 1.0;
    int j2 = 0;
    while ((1.0 + eps_d / 2.0) != 1.0){
        eps_d /= 2.0f;
        j2++;
    }

    std::cout << "SP (32 bit): " << std::endl;
    std::cout << "eps float: " << std::setprecision(50) << eps_f << std::endl;
    std::cout << "Iterations j1 = " << j1 << std::endl << std::endl;

    std::cout << "DP (64 bit): " << std::endl;
    std::cout << "eps double: " << std::setprecision(50) << eps_d << std::endl;
    std::cout << "Iterations j2 = " << j2 << std::endl;

    return 0;
}