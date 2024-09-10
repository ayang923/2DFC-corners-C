#include <stddef.h>
#include <mkl.h>
#include <math.h>
#include <stdio.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

void phi_1D(rd_mat_t x, rd_mat_t *phi_1D_vals);

q_patch_t q_patch_init(M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY, void* phi_param) {
    q_patch_t q_patch;

    q_patch.M_p = M_p;
    q_patch.J = J;
    q_patch.eps_xi_eta = eps_xi_eta;
    q_patch.eps_xy = eps_xy;
    
    q_patch.n_xi = n_xi;
    q_patch.n_eta = n_eta;
    q_patch.xi_start = xi_start;
    q_patch.xi_end = xi_end;
    q_patch.eta_start = eta_start;
    q_patch.eta_end = eta_end;

    q_patch.f_XY = f_XY;

    q_patch.phi_1D = (phi_1D_t) phi_1D;
    q_patch.phi_param = phi_param;

    return q_patch;
}

MKL_INT q_patch_grid_num_el(q_patch_t *q_patch) {
    return (q_patch->n_xi+1)*(q_patch->n_eta+1);
}

void q_patch_evaulate_M_p(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, rd_mat_t *x, rd_mat_t *y) {
    rd_mat_shape(x, xi.rows, xi.columns);
    rd_mat_shape(y, xi.rows, xi.columns);

    q_patch->M_p.M_p_handle(xi, eta, x, y, q_patch->M_p.extra_param);
}
void q_patch_evaulate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t *J_vals) {
    q_patch->J.J_handle(v, J_vals, q_patch->J.extra_param);
}

void q_patch_xi_mesh(q_patch_t *q_patch, rd_mat_t *xi_mesh_vals) {
    rd_mat_shape(xi_mesh_vals, q_patch->n_xi+1, 1);
    rd_linspace(q_patch->xi_start, q_patch->xi_end, q_patch->n_xi+1, xi_mesh_vals);
}

void q_patch_eta_mesh(q_patch_t *q_patch, rd_mat_t *eta_mesh_vals) {
    rd_mat_shape(eta_mesh_vals, q_patch->n_eta+1, 1);
    rd_linspace(q_patch->eta_start, q_patch->eta_end, q_patch->n_eta+1, eta_mesh_vals);
}

void q_patch_xi_eta_mesh(q_patch_t *q_patch, rd_mat_t *XI_vals, rd_mat_t *ETA_vals) {
    double xi_mesh_data[q_patch->n_xi+1];
    rd_mat_t xi_mesh = rd_mat_init_no_shape(xi_mesh_data);
    q_patch_xi_mesh(q_patch, &xi_mesh);

    double eta_mesh_data[q_patch->n_eta+1];
    rd_mat_t eta_mesh = rd_mat_init_no_shape(eta_mesh_data);
    q_patch_eta_mesh(q_patch, &eta_mesh);

    rd_meshgrid(xi_mesh, eta_mesh, XI_vals, ETA_vals);
}

void q_patch_convert_to_XY(q_patch_t *q_patch, rd_mat_t XI, rd_mat_t ETA, rd_mat_t *X_vals, rd_mat_t *Y_vals) {
    q_patch_evaulate_M_p(q_patch, XI, ETA, X_vals, Y_vals);
}

void q_patch_xy_mesh(q_patch_t *q_patch, rd_mat_t *X_vals, rd_mat_t *Y_vals) {
    double XI_data[q_patch_grid_num_el(q_patch)];
    double ETA_data[q_patch_grid_num_el(q_patch)];
    rd_mat_t XI = rd_mat_init_no_shape(XI_data);
    rd_mat_t ETA = rd_mat_init_no_shape(ETA_data);

    q_patch_xi_eta_mesh(q_patch, &XI, &ETA);
    q_patch_convert_to_XY(q_patch, XI, ETA, X_vals, Y_vals);
}

phi_1D_t return_phi_1D(void) {
    return (phi_1D_t) phi_1D;
}

void phi_1D(rd_mat_t x, rd_mat_t *phi_1D_vals) {
    MKL_INT size = x.rows*x.columns;
    double one = 1.0;
    cblas_dcopy(size, &one, 0, phi_1D_vals->mat_data, 1);
    cblas_daxpy(size, -2, x.mat_data, 1, phi_1D_vals->mat_data, 1);
    cblas_dscal(size, 6, phi_1D_vals->mat_data, 1);
    vdErfc(size, phi_1D_vals->mat_data, phi_1D_vals->mat_data);
    cblas_dscal(size, 0.5, phi_1D_vals->mat_data, 1);
}