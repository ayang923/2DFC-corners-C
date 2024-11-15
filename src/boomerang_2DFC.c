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
    return -2.0/3.0*sin(theta*3*M_PI);
}

double l_2(double theta) {
    double alpha = 3.0/2.0;
    double beta = tan(alpha*M_PI/2);
    return beta*sin(theta*2*M_PI);
}

double l_1_prime(double theta) {
    return -2*M_PI*cos(theta*3*M_PI);
}

double l_2_prime(double theta) {
    double alpha = 3.0/2.0;
    double beta = tan(alpha*M_PI/2);
    return 2*M_PI*beta*cos(theta*2*M_PI);
}

double l_1_dprime(double theta) {
    return  6*pow(M_PI, 2)*sin(theta*3*M_PI);
}

double l_2_dprime(double theta) {
    double alpha = 3.0/2.0;
    double beta = tan(alpha*M_PI/2);
    return -4*pow(M_PI, 2)*beta*sin(theta*2*M_PI);
}

int main() {
    MKL_INT d = 7;

    curve_seq_t curve_seq;
    curve_seq_init(&curve_seq);

    curve_t curve_1;
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) l_1, (scalar_func_t) l_2, (scalar_func_t) l_1_prime, (scalar_func_t) l_2_prime, (scalar_func_t) l_1_dprime, (scalar_func_t) l_2_dprime, 100, 1.0/10.0, 1.0/10.0, 0, 0, 0.005);

    c_patch_t c_patches[curve_seq.n_curves];
    s_patch_t s_patches[curve_seq.n_curves];

    rd_mat_t f_arrays[curve_seq_num_f_mats(&curve_seq)];
    double f_points[curve_seq_num_f_mat_points(&curve_seq, d)];

    curve_seq_construct_patches(&curve_seq, s_patches, c_patches, f_arrays, f_points, f, d, 1e-13, 1e-13);

    print_matrix(*(c_patches[0].W.f_XY));
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