#include <stdlib.h>
#include<mkl.h>
#include <math.h>

#include "c_patch_lib.h"
#include "num_linalg_lib.h"
#include "q_patch_lib.h"
#include "fc_lib.h"

void c2_patch_init(c_patch_t *c_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, MKL_INT d, rd_mat_t *f_L, rd_mat_t *f_W) {
    double h_xi = 1.0/(n_xi-1);
    double h_eta = 1.0/(n_eta-1);
    q_patch_init(&(c_patch->W), M_p, J, eps_xi_eta, eps_xy, n_xi, d, 0.0, 1.0, 0.0, (d-1)*h_eta, f_W);
    q_patch_init(&(c_patch->L), M_p, J, eps_xi_eta, eps_xy, d, n_eta-d+1, 0.0, (d-1)*h_xi, (d-1)*h_eta, 1.0, f_L);
    c_patch->c_patch_type = C2;
}

void c1_patch_init(c_patch_t *c_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, MKL_INT d, rd_mat_t *f_L, rd_mat_t *f_W) {
    double h_xi = 1.0/(n_xi-1);
    double h_eta = 1.0/(n_eta-1);
    q_patch_init(&(c_patch->W), M_p, J, eps_xi_eta, eps_xy, (n_xi+1)/2, d, 0.0, 1.0/2.0, 1.0/2.0, 1.0/2.0+(d-1)*h_eta, f_W);
    q_patch_init(&(c_patch->L), M_p, J, eps_xi_eta, eps_xy, d, (n_eta+1)/2+(d-1), 1.0/2.0, 1.0/2.0+(d-1)*h_xi, 0.0, 1.0/2.0+(d-1)*h_eta, f_L);
    c_patch->c_patch_type = C1;
}

void c1_patch_apply_w_W(c_patch_t *c_patch, s_patch_t *window_patch_W) {
    q_patch_apply_w_normalization_xi_left(&(c_patch->W), &(window_patch_W->Q));
}

void c1_patch_apply_w_L(c_patch_t *c_patch, s_patch_t *window_patch_L) {
    q_patch_apply_w_normalization_eta_down(&(c_patch->L), &(window_patch_L->Q));
}

void c2_patch_apply_w_W(c_patch_t *c_patch, s_patch_t *window_patch_W) {
    q_patch_apply_w_normalization_xi_right(&(c_patch->W), &(window_patch_W->Q));
}

void c2_patch_apply_w_L(c_patch_t *c_patch, s_patch_t *window_patch_L) {
    q_patch_apply_w_normalization_eta_up(&(c_patch->L), &(window_patch_L->Q));
}

MKL_INT c2_patch_FC_W_num_el(c_patch_t *c_patch, MKL_INT C, MKL_INT n_r) {
    return (C*n_r+1) * (c_patch->W.n_xi);
}

MKL_INT c2_patch_FC_L_num_el(c_patch_t *c_patch, MKL_INT C, MKL_INT n_r, MKL_INT d) {
    return (C*n_r+1) * (c_patch->L.n_eta+d-1);
}

MKL_INT c_patch_FC_corner_num_el(MKL_INT C, MKL_INT n_r) {
    return pow(C*n_r + 1, 2);
}

MKL_INT c2_patch_FC(c_patch_t *c_patch, MKL_INT C, MKL_INT n_r, MKL_INT d, rd_mat_t A, rd_mat_t Q, q_patch_t* c2_fcont_patch_L, q_patch_t *c2_fcont_patch_W, q_patch_t *c2_fcont_patch_corner, rd_mat_t *f_L, rd_mat_t *f_W, rd_mat_t *f_corner, double *data_stack) {
    f_L->mat_data = data_stack;
    f_W->mat_data = f_L->mat_data + c2_patch_FC_L_num_el(c_patch, C, n_r, d);
    f_corner->mat_data = f_W->mat_data + c2_patch_FC_W_num_el(c_patch, C, n_r);
    
    q_patch_init(c2_fcont_patch_L, c_patch->L.M_p, c_patch->L.J, c_patch->L.eps_xi_eta, c_patch->L.eps_xy, C*n_r+1, c_patch->L.n_eta+d-1, -C*c_patch->L.h_xi, 0, 0, 1, f_L);
    q_patch_init(c2_fcont_patch_W, c_patch->W.M_p, c_patch->W.J, c_patch->W.eps_xi_eta, c_patch->W.eps_xy, c_patch->W.n_xi, C*n_r+1, 0, 1, -C*c_patch->L.h_eta, 0, f_W);
    q_patch_init(c2_fcont_patch_corner, c_patch->W.M_p, c_patch->W.J, c_patch->W.eps_xi_eta, c_patch->W.eps_xy, C*n_r+1, C*n_r+1, -C*c_patch->W.h_xi, 0, -C*c_patch->W.h_eta, 0, f_corner);

    fcont_gram_blend_S(*(c_patch->W.f_XY), d, A, Q, c2_fcont_patch_W->f_XY);

    // computing matrix transpose
    double fcont_W_T_data[c2_fcont_patch_W->f_XY->rows*c2_fcont_patch_W->f_XY->columns];
    mkl_domatcopy('c', 't', c2_fcont_patch_W->f_XY->rows, c2_fcont_patch_W->f_XY->columns, 1.0, c2_fcont_patch_W->f_XY->mat_data, c2_fcont_patch_W->f_XY->rows, fcont_W_T_data, c2_fcont_patch_W->f_XY->columns);

    rd_mat_t fcont_W_T = rd_mat_init(fcont_W_T_data, c2_fcont_patch_W->f_XY->columns, c2_fcont_patch_W->f_XY->rows);

    fcont_gram_blend_S(fcont_W_T, d, A, Q, c2_fcont_patch_corner->f_XY);
    mkl_dimatcopy('c', 't', C*n_r+1, C*n_r+1, 1.0, c2_fcont_patch_corner->f_XY->mat_data, C*n_r+1, C*n_r+1);

    // constructing matrix for L fcont
    double f_L_combined_data[c2_patch_FC_L_num_el(c_patch, C, n_r, d)];
    double *curr_f_W = c_patch->W.f_XY->mat_data;
    double *curr_f_L = c_patch->L.f_XY->mat_data+1;
    double *curr_combined = f_L_combined_data;
    for (int i = 0; i < d; i++) {
        cblas_dcopy(d, curr_f_W, 1, curr_combined, 1);
        curr_combined += d;
        curr_f_W += d;
        cblas_dcopy(c_patch->L.n_eta-1, curr_f_L, 1, curr_combined, 1);
        curr_combined += c_patch->L.n_eta-1;
        curr_f_L += c_patch->L.n_eta;
    }

    mkl_dimatcopy('c', 't', c_patch->L.n_eta+d-1, d, 1.0, f_L_combined_data, c_patch->L.n_eta+d-1, d);
    rd_mat_t f_L_combined = rd_mat_init(f_L_combined_data, d, c_patch->L.n_eta+d-1);
    fcont_gram_blend_S(f_L_combined, d, A, Q, c2_fcont_patch_L->f_XY);

    mkl_dimatcopy('c', 't', C*n_r+1, c_patch->L.n_eta+d-1, 1.0, c2_fcont_patch_L->f_XY->mat_data, C*n_r+1, c_patch->L.n_eta+d-1);
    rd_mat_shape(c2_fcont_patch_L->f_XY, c_patch->L.n_eta+d-1, C*n_r+1);

    return c_patch_FC_corner_num_el(C, n_r) + c2_patch_FC_L_num_el(c_patch, C, n_r, d) + c2_patch_FC_W_num_el(c_patch, C, n_r);
}
