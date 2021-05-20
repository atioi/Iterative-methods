#include <cmath>

double calc_estimator(double *x_next, double *x_before, int size) {

    double estimator = 0.0;
    double max_estimator = 0.0;

    for (int i = 0; i < size; i++) {
        estimator = fabs(x_next[i] - x_before[i]);
        if (estimator > max_estimator)
            max_estimator = estimator;
    }

    return max_estimator;

}
