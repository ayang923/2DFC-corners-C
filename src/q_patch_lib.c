#include <stddef.h>
#include <mkl.h>
#include <math.h>
#include <stdio.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

void phi_1D(rd_mat_t x, rd_mat_t *phi_1D_vals);

void shift_idx_mesh(ri_mat_t *mat, int min_bound, int max_bound) {
    if (mat->mat_data[0] < min_bound) {
        ri_range(min_bound, 1, min_bound+mat->rows-1, mat);
    }
    if (mat->mat_data[mat->rows-1] > max_bound) {
        ri_range(max_bound-mat->rows+1, 1, max_bound, mat);
    }
}

void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY, void* phi_param) {
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

    q_patch->h_xi = (xi_end-xi_start)/n_xi;
    q_patch->h_eta = (eta_end-eta_start)/n_eta;

    q_patch->f_XY = f_XY;

    q_patch->phi_1D = (phi_1D_t) phi_1D;
    q_patch->phi_param = phi_param;
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

void q_patch_evaluate_f(q_patch_t *q_patch, f_handle_t f) {
    double X_data[q_patch_grid_num_el(q_patch)];
    double Y_data[q_patch_grid_num_el(q_patch)];
    rd_mat_t X = rd_mat_init_no_shape(X_data);
    rd_mat_t Y = rd_mat_init_no_shape(Y_data);

    q_patch_xy_mesh(q_patch, &X, &Y);

    f(X, Y, q_patch->f_XY);
}

void q_patch_in_patch(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, ri_mat_t *in_patch_msk) {
    // assumes xi and eta have save shape
    ri_mat_shape(in_patch_msk, xi.rows, xi.columns);

    for (MKL_INT i = 0; i < xi.rows*xi.columns; i++) {
        in_patch_msk->mat_data[i] = xi.mat_data[i] >= q_patch->xi_start && xi.mat_data[i] <= q_patch->xi_end && eta.mat_data[i] >= q_patch->eta_start && eta.mat_data[i] <= q_patch->eta_end;
    }
}

void q_patch_round_boundary_points(q_patch_t *q_patch, rd_mat_t *xi, rd_mat_t *eta) {
    // assumes xi and eta have same shape
    for(MKL_INT i = 0; i < xi->rows*eta->columns; i++) {
        if (fabs(xi->mat_data[i]-q_patch->xi_start) < q_patch->eps_xi_eta) {
            printf("Rounding xi\n");
            xi->mat_data[i] = q_patch->xi_start;
        }
        else if (fabs(xi->mat_data[i]-q_patch->xi_end) < q_patch->eps_xi_eta) {
            xi->mat_data[i] = q_patch->xi_end;
        }
        if (fabs(eta->mat_data[i]-q_patch->eta_start) < q_patch->eps_xi_eta) {
            printf("Rounding eta\n");
            eta->mat_data[i] = q_patch->eta_start;
        }
        else if (fabs(eta->mat_data[i]-q_patch->eta_end) < q_patch->eps_xi_eta) {
            eta->mat_data[i] = q_patch->eta_end;
        }
    }
}

inverse_M_p_return_type_t q_patch_inverse_M_p(q_patch_t *q_patch, double x, double y, rd_mat_t* initial_guesses_xi, rd_mat_t* initial_guesses_eta) {
    rd_mat_t initial_guesses_xi_mat;
    rd_mat_t initial_guesses_eta_mat;
    if (initial_guesses_xi == NULL || initial_guesses_eta == NULL) {
        int N = 20;
        int N_segment = ceil(N/4);

        double initial_guesses_xi_data[N_segment*4];
        double initial_guesses_eta_data[N_segment*4];

        double xi_mesh_data[N_segment+1];
        double eta_mesh_data[N_segment+1];
        rd_mat_t xi_mesh = rd_mat_init(xi_mesh_data, N_segment+1, 1);
        rd_mat_t eta_mesh = rd_mat_init(eta_mesh_data, N_segment+1, 1);

        rd_linspace(q_patch->xi_start, q_patch->xi_end, N_segment+1, &xi_mesh);
        rd_linspace(q_patch->eta_start, q_patch->eta_end, N_segment+1, &eta_mesh);
        
        cblas_dcopy(N_segment, xi_mesh_data, 1, initial_guesses_xi_data, 1);
        cblas_dcopy(N_segment, xi_mesh_data, 1, initial_guesses_xi_data+N_segment, 1);
        cblas_dcopy(N_segment, &(q_patch->xi_start), 0, initial_guesses_xi_data+2*N_segment, 1);
        cblas_dcopy(N_segment, &(q_patch->xi_end), 0, initial_guesses_xi_data+3*N_segment, 1);

        initial_guesses_xi_mat = rd_mat_init(initial_guesses_xi_data, N_segment*4, 1);

        cblas_dcopy(N_segment, &(q_patch->eta_start), 0, initial_guesses_eta_data, 1);
        cblas_dcopy(N_segment, &(q_patch->eta_end), 0, initial_guesses_eta_data+N_segment, 1);
        cblas_dcopy(N_segment, eta_mesh_data, 1, initial_guesses_eta_data+2*N_segment, 1);
        cblas_dcopy(N_segment, eta_mesh_data, 1, initial_guesses_eta_data+3*N_segment, 1);

        initial_guesses_eta_mat = rd_mat_init(initial_guesses_eta_data, N_segment*4, 1);

        initial_guesses_xi = &initial_guesses_xi_mat;
        initial_guesses_eta = &initial_guesses_eta_mat;
    }

    //preallocation for newton's method
    double J_data[4];
    rd_mat_t J_addr = rd_mat_init(J_data, 2, 2);
    double v_prev_data[2];
    rd_mat_t v_prev = rd_mat_init(v_prev_data, 2, 1);
    double v_data[2];
    rd_mat_t v = rd_mat_init(v_data, 2, 1);
    rd_mat_t xi_scalar = rd_mat_init(v_data, 1, 1);
    rd_mat_t eta_scalar = rd_mat_init(v_data+1, 1, 1);

    double f_v_data[4];
    rd_mat_t M_p_v_x = rd_mat_init(f_v_data, 1, 1);
    rd_mat_t M_p_v_y = rd_mat_init(f_v_data+1, 1, 1);

    double v_diff_data[2];
    int ipiv[2];
    double xy_exact_data[2] = {x, y};

    int in_patch;
    ri_mat_t in_patch_scalar = ri_mat_init(&in_patch, 1, 1);

    int converged = 0;

    for (MKL_INT k = 0; k < initial_guesses_xi->columns; k++) {
        v.mat_data[0] = initial_guesses_xi->mat_data[k];
        v.mat_data[1] = initial_guesses_eta->mat_data[k];

        //newton solve with max iterations 100 and error tolerance given by q_patch
        for (MKL_INT i = 0; i < 100; i++) {
            // evaluates difference between solutions in real space
            q_patch_evaulate_M_p(q_patch, rd_mat_init(v.mat_data, 1, 1), rd_mat_init(v.mat_data+1, 1, 1), &M_p_v_x, &M_p_v_y);
            vdSub(2, f_v_data, xy_exact_data, f_v_data);
            
            // convergence threshold
            vdSub(2, v_data, v_prev_data, v_diff_data);
            if (fabs(f_v_data[cblas_idamax(2, f_v_data, 1)]) < q_patch->eps_xy && cblas_dnrm2(2, v_diff_data, 1) < q_patch->eps_xy) {
                converged = 1;
                break;
            }

            v_prev.mat_data[0] = v.mat_data[0];
            v_prev.mat_data[1] = v.mat_data[1];

            //update step
            q_patch_evaulate_J(q_patch, v, &J_addr);
            LAPACKE_dgesv(LAPACK_COL_MAJOR, 2, 1, J_addr.mat_data, 2, ipiv, f_v_data, 2);
            vdSub(2, v_data, f_v_data, v_data);
        }
        
        q_patch_round_boundary_points(q_patch, &xi_scalar, &eta_scalar);
        q_patch_in_patch(q_patch, xi_scalar, eta_scalar, &in_patch_scalar);
        if (converged && in_patch) {
            return (inverse_M_p_return_type_t) {v_data[0], v_data[1], converged};
        }
    }

    return (inverse_M_p_return_type_t) {v_data[0], v_data[1], converged};
}

locally_compute_return_type_t q_patch_locally_compute(q_patch_t *q_patch, double xi, double eta, int M) {
    int in_patch;
    rd_mat_t xi_scalar = rd_mat_init(&xi, 1, 1);
    rd_mat_t eta_scalar = rd_mat_init(&eta, 1, 1);
    ri_mat_t in_patch_scalar = ri_mat_init(&in_patch, 1, 1);
    q_patch_in_patch(q_patch, xi_scalar, eta_scalar, &in_patch_scalar);

    if (!in_patch) {
        return (locally_compute_return_type_t) {NAN, 0};
    }

    int xi_j = (int) ((xi-q_patch->xi_start)/q_patch->h_xi);
    int eta_j = (int) ((eta-q_patch->eta_start)/q_patch->h_eta);

    int half_M = M/2;

    int interpol_xi_j_mesh_data[M];
    int interpol_eta_j_mesh_data[M];
    ri_mat_t interpol_xi_j_mesh = ri_mat_init(interpol_xi_j_mesh_data, M, 1);
    ri_mat_t interpol_eta_j_mesh = ri_mat_init(interpol_eta_j_mesh_data, M, 1);

    if (M%2) {
        ri_range(xi_j-half_M, 1, xi_j+half_M, &interpol_xi_j_mesh);
        ri_range(eta_j-half_M, 1, eta_j+half_M, &interpol_eta_j_mesh);
    }
    else {
        ri_range(xi_j-half_M+1, 1, xi_j+half_M, &interpol_xi_j_mesh);
        ri_range(eta_j-half_M+1, 1, eta_j+half_M, &interpol_eta_j_mesh);
    }

    shift_idx_mesh(&interpol_xi_j_mesh, 0, q_patch->n_xi);
    shift_idx_mesh(&interpol_eta_j_mesh, 0, q_patch->n_eta);

    double interpol_xi_mesh_data[M];
    double interpol_eta_mesh_data[M];
    rd_mat_t interpol_xi_mesh = rd_mat_init(interpol_xi_mesh_data, M, 1);
    rd_mat_t interpol_eta_mesh = rd_mat_init(interpol_eta_mesh_data, M, 1);

    // j*h+start
    for (MKL_INT i = 0; i < M; i++) {
        interpol_xi_mesh_data[i] = interpol_xi_j_mesh_data[i]*q_patch->h_xi + q_patch->xi_start;
        interpol_eta_mesh_data[i] = interpol_eta_j_mesh_data[i]*q_patch->h_eta + q_patch->eta_start;
    }

    double interpol_xi_exact_data[M];
    rd_mat_t interpol_xi_exact = rd_mat_init(interpol_xi_exact_data, M, 1);

    double interpol_val_data[M];
    MKL_INT interpol_xi_idxs[M];

    rd_mat_t interpol_val = rd_mat_init(interpol_val_data, M, 1);
    for (MKL_INT horz_idx = 0; horz_idx < M; horz_idx++) {
        for (MKL_INT i = 0; i < M; i++) {
            interpol_xi_idxs[i] = sub2ind(q_patch->n_eta+1, q_patch->n_xi+1, (sub_t) {interpol_eta_j_mesh_data[horz_idx], interpol_xi_j_mesh_data[i]});
        }
        vdPackV(M, q_patch->f_XY->mat_data, interpol_xi_idxs, interpol_val_data);
        interpol_xi_exact_data[horz_idx] = barylag(interpol_xi_mesh, interpol_val, xi);
    }

    double f_xy = barylag(interpol_eta_mesh, interpol_xi_exact, eta);
    return (locally_compute_return_type_t) {f_xy, 1};
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