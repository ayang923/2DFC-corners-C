#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "curve_seq_lib.h"
#include "s_patch_lib.h"
#include "q_patch_lib.h"
#include "fc2D_lib.h"
#include "fc_lib.h"

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
    double h = 0.01;
    //reading continuation matrices
    MKL_INT d = 4;
    MKL_INT C = 27;
    MKL_INT n_r = 6;

    MKL_INT M = d+3;

    double A_data[fc_A_numel(d, C, n_r)];
    double Q_data[fc_Q_numel(d)];
    rd_mat_t A = rd_mat_init_no_shape(A_data);
    rd_mat_t Q = rd_mat_init_no_shape(Q_data);

    char A_fp[100];
    char Q_fp[100];
    sprintf(A_fp, "fc_data/A_d%d_C%d_r%d.txt", d, C, n_r);
    sprintf(Q_fp, "fc_data/Q_d%d_C%d_r%d.txt", d, C, n_r);
    read_fc_matrix(d, C, n_r, A_fp, Q_fp, &A, &Q);
    curve_seq_t curve_seq;
    curve_seq_init(&curve_seq);

    curve_t curve_1;
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) l_1, (scalar_func_t) l_2, (scalar_func_t) l_1_prime, (scalar_func_t) l_2_prime, (scalar_func_t) l_1_dprime, (scalar_func_t) l_2_dprime, 0, 1.0/10.0, 1.0/10.0, 0, 0, 0.005);

    FC2D(f, h, curve_seq, 1e-14, 1e-14, d, C, n_r, A, Q, M);

    return 0;
}