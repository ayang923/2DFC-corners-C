#ifndef __S_PATCH_LIB_H__
#define __S_PATCH_LIB_H__

#include "q_patch_lib.h"

typedef void(*M_p_general_handle_t) (rd_mat_t xi, rd_mat_t eta, double H, rd_mat_t *x, rd_mat_t *y, void* extra_param);

typedef void (*J_general_handle_t) (rd_mat_t v, double H, rd_mat_t *J_vals, void* extra_param);

typedef struct J_general {
    J_general_handle_t J_general_handle;
    void* extra_param;
} J_general_t;

typedef struct M_p_general {
    M_p_general_handle_t M_p_general_handle;
    void* extra_param;
} M_p_general_t;

typedef struct s_patch {
    q_patch_t q_patch;
    M_p_general_t M_p_general;
    J_general_t J_general;
    double h;
    double H;
} s_patch_t;

void s_patch_init(s_patch_t *s_patch, M_p_general_t M_p_general, J_general_t J_general, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY, double h, void* phi_param);

#endif