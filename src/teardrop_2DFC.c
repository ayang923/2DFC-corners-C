#include <stdio.h>
#include <mkl.h>
#include <math.h>

#include "curve_seq_lib.h"
#include "s_patch_lib.h"
#include "q_patch_lib.h"
#include "fc_lib.h"
#include "r_cartesian_mesh_lib.h"
#include "time.h"

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
    clock_t start, end;
    start = clock();
    double h = 0.005;
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
    curve_seq_add_curve(&curve_seq, &curve_1, (scalar_func_t) l_1, (scalar_func_t) l_2, (scalar_func_t) l_1_prime, (scalar_func_t) l_2_prime, (scalar_func_t) l_1_dprime, (scalar_func_t) l_2_dprime, 0, 1.0/10.0, 1.0/10.0, 0, 0, h);

    c_patch_t c_patches[curve_seq.n_curves];
    s_patch_t s_patches[curve_seq.n_curves];

    rd_mat_t f_arrays[curve_seq_num_f_mats(&curve_seq)];
    double f_points[curve_seq_num_f_mat_points(&curve_seq, d)];

    curve_seq_construct_patches(&curve_seq, s_patches, c_patches, f_arrays, f_points, f, d, 1e-13, 1e-13);

    MKL_INT num_fc_patches = curve_seq_num_FC_mats(&curve_seq);
    q_patch_t fc_patches[num_fc_patches];

    rd_mat_t fc_mats[num_fc_patches];
    double fc_points[curve_seq_num_FC_points(&curve_seq, s_patches, c_patches, C, n_r, d)];


    double x_min = curve_seq.first_curve->l_1(0);
    double x_max = curve_seq.first_curve->l_1(0);
    double y_min = curve_seq.first_curve->l_2(0);
    double y_max = curve_seq.first_curve->l_2(0);

    q_patch_t *curr_q_patch = fc_patches;
    rd_mat_t *curr_fc_mat = fc_mats;
    double *curr_fc_point = fc_points;
    for (int i = 0; i < curve_seq.n_curves; i++) {
        curr_fc_point += s_patch_FC(s_patches+i, C, n_r,d, A, Q, curr_q_patch, curr_fc_mat, curr_fc_point);
        curr_q_patch += 1;
        curr_fc_mat += 1;

        if(c_patches[i].c_patch_type == C2) {
            curr_fc_point += c2_patch_FC(c_patches+i, C, n_r, d, A, Q, curr_q_patch, curr_q_patch+1, curr_q_patch+2, curr_fc_mat, curr_fc_mat+1, curr_fc_mat+2, curr_fc_point);
        }
        curr_q_patch += 3;
        curr_fc_mat += 3;

        for (int j = 1; j <= 4; j++) {
            q_patch_t *prev_q_patch = curr_q_patch-j;
            x_min = MIN(x_min, prev_q_patch->x_min);
            x_max = MAX(x_max, prev_q_patch->x_max);
            y_min = MIN(y_min, prev_q_patch->y_min);
            y_max = MAX(y_max, prev_q_patch->y_max);
        }
    }

    printf("%f, %f, %f, %f\n", x_min, x_max, y_min, y_max);

    double boundary_X_data[curve_seq_boundary_mesh_num_el(&curve_seq, n_r)];
    double boundary_Y_data[curve_seq_boundary_mesh_num_el(&curve_seq, n_r)];
    rd_mat_t boundary_X = rd_mat_init_no_shape(boundary_X_data);
    rd_mat_t boundary_Y = rd_mat_init_no_shape(boundary_Y_data);

    curve_seq_construct_boundary_mesh(&curve_seq, n_r, &boundary_X, &boundary_Y);

    r_cartesian_mesh_obj_t r_cartesian_mesh_obj;

    MKL_INT n_cartesian_mesh = r_cartesian_n_total(x_min-h, x_max+h, y_min-h, y_max+h, h);
    double R_X_data[n_cartesian_mesh];
    double R_Y_data[n_cartesian_mesh];
    MKL_INT in_interior_data[n_cartesian_mesh];
    double f_R_data[n_cartesian_mesh];
    rd_mat_t R_X = rd_mat_init_no_shape(R_X_data);
    rd_mat_t R_Y = rd_mat_init_no_shape(R_Y_data);
    ri_mat_t in_interior = ri_mat_init_no_shape(in_interior_data);
    rd_mat_t f_R = rd_mat_init_no_shape(f_R_data);
    r_cartesian_mesh_init(&r_cartesian_mesh_obj, x_min-h, x_max+h, y_min-h, y_max+h, h, boundary_X, boundary_Y, &R_X, &R_Y, &in_interior, &f_R);

    for (int i = 0; i < num_fc_patches; i++) {
        r_cartesian_mesh_interpolate_patch(&r_cartesian_mesh_obj, fc_patches+i, M);
    }

    r_cartesian_mesh_fill_interior(&r_cartesian_mesh_obj, f);

    FILE *fp;
    fp = freopen("output.txt", "w", stdout);

    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }

    print_matrix(f_R);

    fclose(fp);
    freopen("/dev/tty", "w", stdout);

    end = clock();

    printf("Time: %f\n", ((double) (end-start))/CLOCKS_PER_SEC);
    return 0;
}