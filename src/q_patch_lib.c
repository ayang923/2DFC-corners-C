#include <stddef.h>
#include <mkl.h>
#include <math.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

void phi_1D(rd_mat_t x, rd_mat_t phi_1D_vals);

void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t f_XY, phi_param_t phi_param) {
    q_patch->M_p = M_p;
    q_patch->J = J;
    q_patch->eps_xi_eta = eps_xi_eta;
    q_patch->eps_xy = eps_xy;
    
    q_patch->n_xi = n_xi;
    q_patch->n_eta = n_eta;
    q_patch->xi_start = xi_start;
    q_patch->xi_end = xi_end;
    q_patch->eta_start = eta_start;
    q_patch->eta_end = eta_end;

    q_patch->f_XY = f_XY;

    q_patch->phi_1D = (phi_1D_t) phi_1D;
    q_patch->phi_param = phi_param;
}

void evaulate_M_p(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, rd_mat_t x, rd_mat_t y) {
    q_patch->M_p.M_p_handle(xi, eta, x, y, q_patch->M_p.extra_param);
}
void evaulate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t J_vals) {
    q_patch->J(v, J_vals, NULL);
}

phi_1D_t return_phi_1D(void) {
    return (phi_1D_t) phi_1D;
}

void phi_1D(rd_mat_t x, rd_mat_t phi_1D_vals) {
    MKL_INT size = x.rows*x.columns;
    double one = 1.0;
    cblas_dcopy(size, &one, 0, phi_1D_vals.mat_data, 1);
    cblas_daxpy(size, -2, x.mat_data, 1, phi_1D_vals.mat_data, 1);
    cblas_dscal(size, 6, phi_1D_vals.mat_data, 1);
    vdErfc(size, phi_1D_vals.mat_data, phi_1D_vals.mat_data);
    cblas_dscal(size, 0.5, phi_1D_vals.mat_data, 1);
}