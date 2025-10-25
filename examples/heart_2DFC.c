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

const double ALPH = 1.99;

double f(double x, double y) {
  //    return -(pow(x+1, 2)+pow(y+1, 2))*sin(M_PI*(x-0.1))*cos(M_PI*y);
  //    return -(pow(x+1, 2)+pow(y+1, 2))*sin(M_PI*(x-0.1))*sin(M_PI*y);
  return pow(y, 2)*cos(M_PI*(x-0.1));
}

double l_1(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return bet*cos((1+ALPH)*M_PI*theta)-sin((1+ALPH)*M_PI*theta)-bet;
}

double l_2(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return bet*sin((1+ALPH)*M_PI*theta)+cos((1+ALPH)*M_PI*theta)-cos(M_PI*theta);
}

double l_1_prime(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return -bet*(1+ALPH)*M_PI*sin((1+ALPH)*M_PI*theta)-(1+ALPH)*M_PI*cos((1+ALPH)*M_PI*theta);
}

double l_2_prime(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return bet*(1+ALPH)*M_PI*cos((1+ALPH)*M_PI*theta)-(1+ALPH)*M_PI*sin((1+ALPH)*M_PI*theta)+M_PI*sin(M_PI*theta);
}

double l_1_dprime(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return -bet*pow((1+ALPH)*M_PI, 2)*cos((1+ALPH)*M_PI*theta)+pow((1+ALPH)*M_PI, 2)*sin((1+ALPH)*M_PI*theta);
}

double l_2_dprime(double theta) {
    double bet = tan(ALPH*M_PI/2);
    return -bet*pow((1+ALPH)*M_PI, 2)*sin((1+ALPH)*M_PI*theta)-pow((1+ALPH)*M_PI, 2)*cos((1+ALPH)*M_PI*theta) + pow(M_PI, 2)*cos(M_PI*theta);

}

int main() {
    //reading continuation matrices
    MKL_INT d = 7;
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


    double n_frac_c = 0.1;
    double n_frac_S = 0.8;

    double h = 0.0005;
    double h_tan = 2.1*h;
    double h_norm = h_tan;
    MKL_INT n_curve = ceil(3.5/h_norm);
    
    curve_seq_t curve_seq;
    curve_seq_init(&curve_seq);
    
    curve_t curve_1;
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) l_1, (scalar_func_t) l_2, (scalar_func_t) l_1_prime, (scalar_func_t) l_2_prime, (scalar_func_t) l_1_dprime, (scalar_func_t) l_2_dprime, n_curve, n_frac_c, n_frac_c, n_frac_S, n_frac_S, h_tan);

    FC2D(f, h, curve_seq, 1e-14, 1e-14, d, C, n_r, A, Q, M);

    return 0;
}
