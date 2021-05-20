#include "other.h"
#include "estimator.h"
#include "residuum.h"
#include <iostream>

void Jacobi(double *A, int size, double *b, double *x_initial) {

    double M[size * size];
    double C[size];

    /* Wzór operacyjny:
     *
     *  x_next = - 1 / D * ( L + U ) * x_before + 1 / D * b
     *  M = - 1 / D * (L+U)
     *  C = 1/D * b
     *  x_next = M * x_before + C
     *
     * */

    for (int n = 0; n < size; n++)
        C[n] = b[n] / A[n * (size + 1)];

    for (int k = 0; k < size; k++)
        for (int j = 0; j < size; j++) {
            int diag_index = k * (size + 1);
            int i = j + (k * size);
            M[i] = -A[i] / A[diag_index];
            M[diag_index] = 0;
        }

    std::cout << std::endl << "M matrix:" << std::endl;
    print_matrix(M, size);
    std::cout << std::endl << "C vector: " << std::endl;
    print_vector(C, size);

    // Parameters:
    int n_max = 10; // Maximum number of iterations
    double TOLX = 1e-9; // Given error tolerance
    double TOLF = 1e-9; // Given residuum tolerance

    double residuum = 0.0;
    double estimator = 0.0;

    double x_next[size];
    double *x_before = x_initial;


    // Main loop:
    for (int n = 0; n < n_max; n++) {


        // Calculate x_next:
        for (int k = 0; k < size; k++) {
            x_next[k] = 0;
            for (int j = 0; j < size; j++) {
                int i = j + (k * size);
                x_next[k] += (M[i] * x_before[j]);
            }
            x_next[k] += C[k];
        }

        residuum = calc_residuum(A, b, x_next, size);
        estimator = calc_estimator(x_next, x_before, size);
        print_results(n, x_next, size, residuum, estimator);

        if (residuum <= TOLF && estimator <= TOLX)
            break;

        // Copy x_next to x_before:
        for (int i = 0; i < size; i++)
            x_before[i] = x_next[i];

    }
}
