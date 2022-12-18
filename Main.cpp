#define WITH_NEON
#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include "Matrix.hpp"

int main(){
    try{
        Matrix<float> a(4,4);
        Matrix<float> b(4,4);
        fillMatrix(a);
        fillMatrix(b);
        printMatrix(a*b);
    }
    catch(...){}
}
