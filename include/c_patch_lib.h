#ifndef __C_PATCH_LIB_H__
#define __C_PATCH_LIB_H__

#include <mkl.h>
#include <stdlib.h>

#include "num_linalg_lib.h"
#include "q_patch_lib.h"
#include "s_patch_lib.h"

typedef enum c_patch_type {
    C1,
    C2,
} c_patch_type_t;

typedef struct c_patch {
    q_patch_t L;
    q_patch_t W;
    c_patch_type_t c_patch_type;
} c_patch_t;

// void c1_patch_init(c_patch_t *c_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, MKL_INT d, rd_mat_t *f_L, rd_mat_t *f_W);

void c1_patch_init(c_patch_t *c_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, MKL_INT d, rd_mat_t *f_L, rd_mat_t *f_W);

void c2_patch_init(c_patch_t *c_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, MKL_INT d, rd_mat_t *f_L, rd_mat_t *f_W);

void c1_patch_apply_w_W(c_patch_t *c_patch, s_patch_t *window_patch_W);

void c1_patch_apply_w_L(c_patch_t *c_patch, s_patch_t *window_patch_L);

void c2_patch_apply_w_W(c_patch_t *c_patch, s_patch_t *window_patch_W);

void c2_patch_apply_w_L(c_patch_t *c_patch, s_patch_t *window_patch_L);
#endif