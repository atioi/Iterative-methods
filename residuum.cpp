#include <cmath>


double calc_residuum(double *A, double *b, double *x, int size) {

    double residuum = 0.0;
    double max_residuum = 0.0;

    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++)
            residuum += x[k] * A[i * size + k];

        residuum = b[i] - residuum;

        if (residuum > max_residuum)
            max_residuum = residuum;

        residuum = 0;
    }

    return max_residuum;

}
