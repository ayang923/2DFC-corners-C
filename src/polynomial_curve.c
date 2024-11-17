#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "curve_seq_lib.h"
#include "s_patch_lib.h"
#include "q_patch_lib.h"
#include "fc_lib.h"
#include "r_cartesian_mesh_lib.h"
#include "time.h"
#include "fc2D_lib.h"

double f(double x, double y) {return 4+(1+pow(x, 2) + pow(y, 2))*(sin(2.5*M_PI*x-0.5)+cos(2*M_PI*y-0.5));}

double curve_1_l_1(double theta) {return -theta+1;}

double curve_1_l_2(double theta) {return -theta+1;}

double curve_1_l_1_prime(double theta) {return -1;}

double curve_1_l_2_prime(double theta) {return -1;}

double curve_1_l_1_dprime(double theta) {return 0;}

double curve_1_l_2_dprime(double theta) {return 0;}

double curve_2_l_1(double theta) {return theta;}

double curve_2_l_2(double theta) {return pow(theta, 3);}

double curve_2_l_1_prime(double theta) {return 1;}

double curve_2_l_2_prime(double theta) {return 3*pow(theta, 2);}

double curve_2_l_1_dprime(double theta) {return 0;}

double curve_2_l_2_dprime(double theta) {return 6*theta;}



int main() {

    double h = 0.0005;
    double h_tan = 0.001;
    //reading continuation matrices
    MKL_INT d = 4;
    MKL_INT C = 27;
    MKL_INT n_r = 6;

    MKL_INT M = d+3;

    double A_data[fc_A_numel(d, C, n_r)];
    double Q_data[fc_Q_numel(d)];
    rd_mat_t A = rd_mat_init_no_shape(A_data);
    rd_mat_t Q = rd_mat_init_no_shape(Q_data);

    read_fc_matrix(d, C, n_r, "fc_data/A_d4_C27_r6.txt", "fc_data/Q_d4_C27_r6.txt", &A, &Q);

    curve_seq_t curve_seq;
    curve_seq_init(&curve_seq);

    curve_t curve_1;
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) curve_1_l_1, (scalar_func_t) curve_1_l_2, (scalar_func_t) curve_1_l_1_prime, (scalar_func_t) curve_1_l_2_prime, (scalar_func_t) curve_1_l_1_dprime, (scalar_func_t) curve_1_l_2_dprime, 0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, h_tan);

    curve_t curve_2;
    curve_seq_add_curve(&curve_seq, &curve_2, (scalar_func_t) curve_2_l_1, (scalar_func_t) curve_2_l_2, (scalar_func_t) curve_2_l_1_prime, (scalar_func_t) curve_2_l_2_prime, (scalar_func_t) curve_2_l_1_dprime, (scalar_func_t) curve_2_l_2_dprime, 0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, h_tan);

    FC2D(f, h, curve_seq, 1e-13, 1e-13, d, C, n_r, A, Q, M);

    return 0;
}