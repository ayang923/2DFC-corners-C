#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "q_patch_lib.h"
#include "s_patch_lib.h"
#include "num_linalg_lib.h"
#include "fc_lib.h"

typedef struct M_p_S_extra_param {
    double theta_A;
    double theta_B;
} M_p_S_extra_param_t;

typedef struct J_S_extra_param {
    double theta_A;
    double theta_B;
} J_S_extra_param_t;

void M_p_general_S(rd_mat_t xi, rd_mat_t eta, double H, rd_mat_t *x, rd_mat_t *y, void* extra_param) {
    M_p_S_extra_param_t* S_extra_param = (M_p_S_extra_param_t*) extra_param;

    // assumes xi and eta are of same size, using for loopos for ease of implementation
    for (MKL_INT i = 0; i < xi.rows*xi.columns; i++) {
        double l_A = xi.mat_data[i]*(S_extra_param->theta_B - S_extra_param->theta_A) + S_extra_param->theta_A;
        double nu_1 = -0.7*cos(l_A) / sqrt(pow(0.7*cos(l_A), 2) + pow(sin(l_A)+0.7*sin(2*l_A), 2));
        double nu_2 = (-sin(l_A)-0.7*sin(2*l_A)) / sqrt(pow(0.7*cos(l_A), 2) + pow(sin(l_A)+0.7*sin(2*l_A), 2));

        x->mat_data[i] = cos(l_A) + 0.35*cos(2*l_A)-0.35 + eta.mat_data[i] * H * nu_1;
        y->mat_data[i] = 0.7*sin(l_A) + eta.mat_data[i] * H * nu_2;
    }
}

void f(rd_mat_t x, rd_mat_t y, rd_mat_t *f_xy) {
    rd_mat_shape(f_xy, x.rows, x.columns);
    for (MKL_INT i = 0; i < x.rows*x.columns; i++) {
        f_xy->mat_data[i] = -(pow(x.mat_data[i], 6) + pow(y.mat_data[i], 6))*sin(10*M_PI*x.mat_data[i])*sin(10*M_PI*y.mat_data[i]);
    }
}

int main() {
    M_p_S_extra_param_t M_p_S1_extra_param = {0, M_PI};
    M_p_S_extra_param_t M_p_S2_extra_param = {M_PI, 2*M_PI};
    M_p_general_t M_p_general_S1 = {(M_p_general_handle_t) M_p_general_S, (void *) &M_p_S1_extra_param};
    M_p_general_t M_p_general_S2 = {(M_p_general_handle_t) M_p_general_S, (void *) &M_p_S2_extra_param};

    double h = 0.01;
    MKL_INT n_S = (MKL_INT) round(M_PI/h);

    MKL_INT d = 4;
    MKL_INT C = 27;
    MKL_INT n_r = 6;

    double f_S1_data[(n_S+1)*(d+1)];
    rd_mat_t f_S1 = rd_mat_init(f_S1_data, d+1, n_S+1);
    double f_S2_data[(n_S+1)*(d+1)];
    rd_mat_t f_S2 = rd_mat_init(f_S2_data, d+1, n_S+1);

    s_patch_t s1_patch;
    s_patch_init(&s1_patch, M_p_general_S1, (J_general_t) {(J_general_handle_t) M_p_general_S, NULL}, 1e-13, 1e-13, n_S, d, 0, 1, 0, 1, &f_S1, h, NULL);
    s_patch_t s2_patch;
    s_patch_init(&s2_patch, M_p_general_S2, (J_general_t) {(J_general_handle_t) M_p_general_S, NULL}, 1e-13, 1e-13, n_S, d, 0, 1, 0, 1, &f_S2, h, NULL);

    q_patch_evaluate_f((q_patch_t*) &s1_patch, (f_handle_t) f);
    q_patch_evaluate_f((q_patch_t*) &s2_patch, (f_handle_t) f);
    
    //reading continuation matrices
    double A_data[fc_A_numel(d, C, n_r)];
    double Q_data[fc_Q_numel(d)];
    rd_mat_t A = {A_data, 0, 0};
    rd_mat_t Q = {Q_data, 0, 0};

    read_fc_matrix(d, C, n_r, "fc_data/A_d4_C27_r6.txt", "fc_data/Q_d4_C27_r6.txt", &A, &Q);

    s_patch_t s1_FC_patch;
    double s1_f_FC_data[s_patch_FC_num_el(&s1_patch, C, n_r)];
    rd_mat_t s1_f_FC = {s1_f_FC_data, 0, 0};

    s_patch_t s2_FC_patch;
    double s2_f_FC_data[s_patch_FC_num_el(&s2_patch, C, n_r)];
    rd_mat_t s2_f_FC = {s2_f_FC_data, 0, 0};

    s_patch_FC_init(&s1_patch, C, n_r, &s1_f_FC, &s1_FC_patch);
    s_patch_FC_init(&s2_patch, C, n_r, &s2_f_FC, &s2_FC_patch);

    s_patch_FC(&s1_patch, d, A, Q, NULL, &s1_FC_patch);
    s_patch_FC(&s2_patch, d, A, Q, NULL, &s2_FC_patch);

    return 0;
}