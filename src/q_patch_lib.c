#include <stddef.h>
#include <cblas.h>
#include <math.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

void phi_1D(rd_mat_t x, rd_mat_t phi_1D_vals);

void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, size_t n_xi, size_t n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t f_XY, phi_param_t phi_param) {
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
    q_patch->M_p(xi, eta, x, y, NULL);
}
void evaulate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t J_vals) {
    q_patch->J(v, J_vals, NULL);
}

phi_1D_t return_phi_1D(void) {
    return (phi_1D_t) phi_1D;
}

void phi_1D(rd_mat_t x, rd_mat_t phi_1D_vals) {
    size_t size = x.rows*x.columns;
    cblas_dcopy(size, x.mat_data, 1, phi_1D_vals.mat_data, 1);
    cblas_dscal(size, -2, phi_1D_vals.mat_data, 1);

    double mat_data[size];
    rd_mat_t ones_mat = {mat_data, x.rows, x.columns};
    rd_ones(ones_mat);

    cblas_daxpy(size, 1, ones_mat.mat_data, 1, phi_1D_vals.mat_data, 1);
    cblas_dscal(size, 6, phi_1D_vals.mat_data, 1);

    for (size_t i = 0; i < size; i++) {
        phi_1D_vals.mat_data[i] = erfc(phi_1D_vals.mat_data[i]);
    }

    cblas_dscal(size, 0.5, phi_1D_vals.mat_data, 1);
}