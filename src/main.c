#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

typedef struct M_p_C2_extra_param {
    double theta_A;
    double theta_B;
    double theta_C;
} M_p_C2_extra_param_t;

void M_p_C2(rd_mat_t xi, rd_mat_t eta, rd_mat_t x, rd_mat_t y, void* extra_param) {
    M_p_C2_extra_param_t* C2_extra_param = (M_p_C2_extra_param_t*) extra_param;

    // assumes xi and eta are column vectors of same side, using for loopos for ease of implementation
    for (MKL_INT i = 0; i < xi.rows; i++) {
        double l_A = xi.mat_data[i]*(C2_extra_param->theta_A - C2_extra_param->theta_C) + C2_extra_param->theta_C;
        double l_B = eta.mat_data[i]*(C2_extra_param->theta_B - C2_extra_param->theta_C - 2*M_PI) + C2_extra_param->theta_C + 2*M_PI;

        // printf("(%f, %f)\n", l_A, l_B);

        x.mat_data[i] = 2*sin(l_A/2) + 2*sin(l_B/2) - 2*sin(C2_extra_param->theta_C/2);
        y.mat_data[i] = -sin(l_A) - sin(l_B) - sin(C2_extra_param->theta_C);
    }
}

int main() {

    double xi_mesh_data[21];
    rd_mat_t xi_mesh = {xi_mesh_data, 21, 1};
    rd_linspace(0, 1, 21, xi_mesh);

    double eta_mesh_data[41];
    rd_mat_t eta_mesh = {eta_mesh_data, 41, 1};
    rd_linspace(0, 1, 41, eta_mesh);
    // for (MKL_INT i = 0; i < xi_mesh.rows; i++) {
    //     printf("%f\n", xi_mesh.mat_data[i]);
    // }

    double XI_data[xi_mesh.rows * eta_mesh.rows];
    double ETA_data[xi_mesh.rows * eta_mesh.rows];
    rd_mat_t XI = {XI_data, 0, 0};
    rd_mat_t ETA = {ETA_data, 0, 0};

    rd_meshgrid(xi_mesh, eta_mesh, &XI, &ETA);
    printf("%d, %d\n", XI.rows, XI.columns);
    for (int i = 0; i < XI.rows; i++) {
        for (int j = 0; j < XI.columns; j++) {
            MKL_INT idx = sub2ind(XI.rows, XI.columns, (sub_t) {i, j});
            printf("%f ", ETA.mat_data[idx]);
        }
        printf("\n");
    }

    M_p_C2_extra_param_t C2_params = {0.4, 2*M_PI-0.4, 0};
    double x[21];
    rd_mat_t x_mesh = {x, 21, 1};
    double y[21];
    rd_mat_t y_mesh = {y, 21, 1};

    M_p_C2(xi_mesh, xi_mesh, x_mesh, y_mesh, (void*) &C2_params);

    // for (MKL_INT i = 0; i < 21; i++) {
    //     // printf("(%f, %f)\n", xi_mesh.mat_data[i], xi_mesh.mat_data[i]);
    //     printf("(%f, %f)\n", x[i], y[i]);
    // }

    return 0;
}