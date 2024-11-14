#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "curve_seq_lib.h"
#include "s_patch_lib.h"
#include "q_patch_lib.h"

double f(double x, double y) {
    return 4+(1+pow(x, 2) + pow(y, 2))*(sin(2.5*M_PI*x-0.5)+cos(2*M_PI*y-0.5));
}

double l_1(double theta) {
    return 2*sin(theta*M_PI);
}

double l_2(double theta) {
    return -sin(theta*2*M_PI);
}

double l_1_prime(double theta) {
    return 2*M_PI*cos(theta*M_PI);
}

double l_2_prime(double theta) {
    return -2*M_PI*cos(theta*2*M_PI);
}

double l_1_dprime(double theta) {
    return -2*pow(M_PI, 2)*sin(theta*M_PI);
}

double l_2_dprime(double theta) {
    return 4*pow(M_PI, 2)*sin(theta*2*M_PI);
}

int main() {
    MKL_INT d = 7;

    curve_seq_t curve_seq;
    curve_seq_init(&curve_seq);

    curve_t curve_1;
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) l_1, (scalar_func_t) l_2, (scalar_func_t) l_1_prime, (scalar_func_t) l_2_prime, (scalar_func_t) l_1_dprime, (scalar_func_t) l_2_dprime, 30, 1.0/10.0, 1.0/10.0, 0, 0, 0.005);

    s_patch_t s_patch;
    double s_patch_f_XY_data[curve_S_patch_num_el(&curve_1, 7)];
    rd_mat_t s_patch_f_XY = rd_mat_init_no_shape(s_patch_f_XY_data);
    curve_construct_S_patch(&curve_1, &s_patch, &s_patch_f_XY, (scalar_func_2D_t) f, d, 1e-13, 1e-13);

    inverse_M_p_return_type_t result = q_patch_inverse_M_p(&(s_patch.q_patch), 0.65, 0.60, NULL, NULL);
    printf("%lf, %lf, %d\n", result.xi, result.eta, result.converged);

    // M_p_S_extra_param_t M_p_S1_extra_param = {0, M_PI};
    // M_p_S_extra_param_t M_p_S2_extra_param = {M_PI, 2*M_PI};
    // M_p_general_t M_p_general_S1 = {(M_p_general_handle_t) M_p_general_S, (void *) &M_p_S1_extra_param};
    // M_p_general_t M_p_general_S2 = {(M_p_general_handle_t) M_p_general_S, (void *) &M_p_S2_extra_param};

    // double h = 0.01;
    // MKL_INT n_S = (MKL_INT) round(M_PI/h);

    // MKL_INT d = 4;
    // MKL_INT C = 27;
    // MKL_INT n_r = 6;

    // double f_S1_data[(n_S+1)*(d+1)];
    // rd_mat_t f_S1 = rd_mat_init(f_S1_data, d+1, n_S+1);
    // double f_S2_data[(n_S+1)*(d+1)];
    // rd_mat_t f_S2 = rd_mat_init(f_S2_data, d+1, n_S+1);

    // s_patch_t s1_patch;
    // s_patch_init(&s1_patch, M_p_general_S1, (J_general_t) {(J_general_handle_t) M_p_general_S, NULL}, 1e-13, 1e-13, n_S, d, 0, 1, 0, 1, &f_S1, h, NULL);
    // s_patch_t s2_patch;
    // s_patch_init(&s2_patch, M_p_general_S2, (J_general_t) {(J_general_handle_t) M_p_general_S, NULL}, 1e-13, 1e-13, n_S, d, 0, 1, 0, 1, &f_S2, h, NULL);

    // q_patch_evaluate_f((q_patch_t*) &s1_patch, (f_handle_t) f);
    // q_patch_evaluate_f((q_patch_t*) &s2_patch, (f_handle_t) f);
    
    // //reading continuation matrices
    // double A_data[fc_A_numel(d, C, n_r)];
    // double Q_data[fc_Q_numel(d)];
    // rd_mat_t A = {A_data, 0, 0};
    // rd_mat_t Q = {Q_data, 0, 0};

    // read_fc_matrix(d, C, n_r, "fc_data/A_d4_C27_r6.txt", "fc_data/Q_d4_C27_r6.txt", &A, &Q);

    // s_patch_t s1_FC_patch;
    // double s1_f_FC_data[s_patch_FC_num_el(&s1_patch, C, n_r)];
    // rd_mat_t s1_f_FC = {s1_f_FC_data, 0, 0};

    // s_patch_t s2_FC_patch;
    // double s2_f_FC_data[s_patch_FC_num_el(&s2_patch, C, n_r)];
    // rd_mat_t s2_f_FC = {s2_f_FC_data, 0, 0};

    // s_patch_FC_init(&s1_patch, C, n_r, &s1_f_FC, &s1_FC_patch);
    // s_patch_FC_init(&s2_patch, C, n_r, &s2_f_FC, &s2_FC_patch);

    // s_patch_FC(&s1_patch, d, A, Q, NULL, &s1_FC_patch);
    // s_patch_FC(&s2_patch, d, A, Q, NULL, &s2_FC_patch);

    return 0;
}