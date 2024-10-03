
#include "s_patch_lib.h"
#include  "q_patch_lib.h"

void s_M_p_wrapper_function(rd_mat_t xi, rd_mat_t eta, rd_mat_t *x, rd_mat_t *y, void* extra_param) {
    s_patch_t *s_patch = (s_patch_t*) extra_param;
    s_patch->M_p_general.M_p_general_handle(xi, eta, s_patch->H, x, y, s_patch->M_p_general.extra_param);
}

void s_J_wrapper_function(rd_mat_t v, rd_mat_t *J_vals, void* extra_param) {
    s_patch_t *s_patch = (s_patch_t*) extra_param;
    s_patch->J_general.J_general_handle(v, s_patch->H, J_vals, s_patch->J_general.extra_param);
}

void s_patch_init(s_patch_t *s_patch, M_p_general_t M_p_general, J_general_t J_general, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY, double h, void* phi_param) {
    q_patch_init(&(s_patch->q_patch), (M_p_t) {(M_p_handle_t) s_M_p_wrapper_function, (void*) s_patch}, (J_t) {(J_handle_t) s_J_wrapper_function, (void*) s_patch}, eps_xi_eta, eps_xy, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi_param);
    
    s_patch->h = h;
    s_patch->H = h / s_patch->q_patch.h_eta;
    s_patch->M_p_general = M_p_general;
    s_patch->J_general = J_general;
}

// void s_patch_FC(MKL_INT C, MKL_INT n_r, MKL_INT d, rd_mat)