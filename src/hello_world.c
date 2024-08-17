#include <stdio.h>
#include "q_patch_lib.h"
#include <cblas.h>

int main() {
    // Define a 3x3 matrix and a 3x1 vector
    double A[9] = {1.0, 2.0, 3.0,
                   4.0, 5.0, 6.0,
                   7.0, 8.0, 9.0};
    double x[3] = {1.0, 0.0, -1.0};
    double y[3] = {0.0, 0.0, 0.0};  // Result vector

    // Perform the matrix-vector multiplication y = A * x
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, A, 3, x, 1, 0.0, y, 1);

    // Print the result
    printf("Resulting vector y:\n");
    for (int i = 0; i < 3; i++) {
        printf("%f\n", y[i]);
    }

    return 0;
}