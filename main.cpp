#include <iostream>
#include "methods.h"

int main() {
    double A[] = {
            100.0, -1.0, 2.0, -3.0,
            1.0, 200.0, -4.0, 5.0,
            -2.0, 4.0, 300.0, -6.0,
            3.0, -5.0, 6.0, 400.0
    };
    double b[] = {116.0, -226.0, 912.0, -1174.0};
    double x_initial[] = {2.0, 2.0, 2.0, 2.0};

    double relaxation_factor = 0.5;
    int size = 4;

//    Jacobi(A, size, b, x_initial);
//    Gauss_Seidel(A, size, b, x_initial);
    SOR(A, size, b, relaxation_factor, x_initial);


}
