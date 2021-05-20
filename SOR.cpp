#include "other.h"
#include <iostream>
#include "estimator.h"
#include "residuum.h"


void SOR(double *A, int size, double *b, double relaxation_factor, double *x_initial) {


    /*
     * Wzór operacyjny:
     *
     * ( L + 1 / relaxation_factor * D) * x_next = - [ ( 1 - 1 / relaxation_factor ) * D + U ] * x_before + b
     *  ^---------  LEFT ------------^            ^--------------------- RIGHT ---------------^
     *
     * */

    double RIGHT[size * size];
    double LEFT[size * size];

    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            int index = row * size + col;

            LEFT[index] = 0;
            RIGHT[index] = 0;

            if (row > col) {
                // Creates LEFT:
                LEFT[index] = A[index];

            } else if (row < col) {
                // Creates RIGHT:
                RIGHT[index] = -A[index];

            } else {
                LEFT[index] = 1 / relaxation_factor * A[index];
                RIGHT[index] = -(1 - 1 / relaxation_factor) * A[index];
            }
        }
    }

    std::cout << std::endl << "RIGHT:" << std::endl;
    print_matrix(RIGHT, size);
    std::cout << std::endl << "LEFT:" << std::endl;
    print_matrix(LEFT, size);


    // Parameters:
    int n_max = 40; // Maximum number of iterations
    double TOLX = 1e-9; // Given error tolerance
    double TOLF = 1e-9; // Given residuum tolerance

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
                sum += x_before[c] * RIGHT[r * size + c];
            }
            temp[r] = b[r] + sum;
        }

        x_next[0] = temp[0] / LEFT[0];

        for (int r = 1; r < size; r++) {
            sum = 0.0;
            for (int c = 0; c < r; c++) {
                sum += LEFT[r * size + c] * x_next[c];
            }
            x_next[r] = (temp[r] - sum) / LEFT[r * (size + 1)];
        }

        std::cout << std::endl;
        std::cout << std::endl;
        print_vector(x_next, size);
        std::cout << std::endl;


        residuum = calc_residuum(A, b, x_next, size);
        estimator = calc_estimator(x_next, x_before, size);
        print_results(n, x_next, size, residuum, estimator);

        if (residuum <= TOLF && estimator <= TOLX)
            break;

        for (int i = 0; i < size; i++)
            x_before[i] = x_next[i];


    }


}
