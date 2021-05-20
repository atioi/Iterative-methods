#include <iostream>
#include <iomanip>

void print_matrix(double M[], int size) {
    for (int r = 0; r < size; r++) {
        for (int c = 0; c < size; c++) {
            int index = r * size + c;
            std::cout << std::setw(10) << M[index];
        }
        std::cout << std::endl;
    }
}

void print_vector(double V[], int size) {
    for (int i = 0; i < size; i++)
        std::cout << V[i] << ' ';
}

void print_results(int iteration, double x_next[], int size, double residuum, double estimator) {
    std::cout << std::endl << std::endl
              << "Iteration: "
              << iteration
              << std::endl;
    std::cout << '[';
    print_vector(x_next, size);
    std::cout << ']';
    std::cout << "  Residuum: " << residuum << "  Estimator: " << estimator;
}
