#include <stdlib.h>
#include<mkl.h>

#include "c_patch_lib.h"
#include "num_linalg_lib.h"
#include "q_patch_lib.h"

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