#include<iostream>
#include<cmath>

double power(double base, int exponent){
    double result = 1;
    for(int i = 0; i<exponent; i++){
        result = result*base;
    }
    return result;
}

void print_pow(double base, int exponent){ //void = function with no return
    double mypower = power(base, exponent);
    std::cout << base << " raised to the power of " << exponent<< " is " << mypower << std::endl;
}

int main(){
    double base;
    int exponent;
    std::cout << "What is the base?" << std::endl;
    std::cin >> base;
    std::cout << "What is the exponenet?" << std::endl;
    std::cin >> exponent;
    print_pow(base, exponent);
}

