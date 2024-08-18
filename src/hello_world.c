#include <stdio.h>
#include <cblas.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

int main() {
    phi_1D_t phi_1D = return_phi_1D();

    double mat_data_x[2] = {0.5, 0.25};
    rd_mat_t x = {mat_data_x, 2, 1};
    double mat_data_phi[2];
    rd_mat_t phi_1D_vals = {mat_data_phi, 2, 1};

    phi_1D(x, phi_1D_vals);

    for (size_t i = 0; i < phi_1D_vals.rows; i++) {
        printf("%f\n", phi_1D_vals.mat_data[i]);
    }
    return 0;
}