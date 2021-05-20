#include "other.h"
#include "residuum.h"
#include "estimator.h"
#include <iostream>


void Gauss_Seidel(double *A, int size, double *b, double *x_initial) {

    /*
     *
     * Wzór operacyjny:
     * ( L + D ) * n_next = - U * x_before + b
     *
     * */

    double LD[size * size]; // = L+D
    double U[size * size]; // = U


    // Creates LD and U matrix
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            int index = row * size + col;

            if (row >= col) {
                // LD matrix:
                LD[index] = A[index];
                U[index] = 0;
            } else {
                // U matrix:
                U[index] = A[index];
                LD[index] = 0;
            }
        }
    }


    std::cout << std::endl << "L + D matrix:" << std::endl;
    print_matrix(LD, size);
    std::cout << std::endl << "U matrix:" << std::endl;
    print_matrix(U, size);

    // Parameters:
    int n_max = 10; // Maximum number of iterations
    double TOLX = 1e-3; // Given error tolerance
    double TOLF = 1e-3; // Given residuum tolerance

    double residuum = 0.0;
    double estimator = 0.0;

    double x_next[size];
    double *x_before = x_initial;
    double temp[size];
    double sum = 0.0;

    for (int n = 0; n < n_max; n++) {

        for (int r = 0; r < size; r++) {
            sum = 0.0;
            for (int c = 0; c < size; c++) {
                sum += x_before[c] * U[r * size + c];
            }
            temp[r] = b[r] - sum;
        }

        x_next[0] = temp[0] / LD[0];

        for (int r = 1; r < size; r++) {
            sum = 0.0;
            for (int c = 0; c < r; c++) {
                sum += LD[r * size + c] * x_next[c];
            }
            x_next[r] = (temp[r] - sum) / LD[r * (size + 1)];
        }



        residuum = calc_residuum(A, b, x_next, size);
        estimator = calc_estimator(x_next, x_before, size);
        print_results(n, x_next, size, residuum, estimator);

        if (residuum <= TOLF && estimator <= TOLX)
            break;

        for (int i = 0; i < size; i++)
            x_before[i] = x_next[i];


    }
};