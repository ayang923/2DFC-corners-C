#include <stddef.h>
#include <mkl.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "q_patch_lib.h"
#include "num_linalg_lib.h"

void w_1D(rd_mat_t x, rd_mat_t *w_1D_vals);

void shift_idx_mesh(ri_mat_t *mat, int min_bound, int max_bound) {
    if (mat->mat_data[0] < min_bound) {
        ri_range(min_bound, 1, min_bound+mat->rows-1, mat);
    }
    if (mat->mat_data[mat->rows-1] > max_bound) {
        ri_range(max_bound-mat->rows+1, 1, max_bound, mat);
    }
}

void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY) {
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

    q_patch->h_xi = (xi_end-xi_start)/(n_xi-1);
    q_patch->h_eta = (eta_end-eta_start)/(n_eta-1);

    q_patch->f_XY = f_XY;
    rd_mat_shape(f_XY, n_eta, n_xi);

    q_patch->w_1D = (w_1D_t) w_1D;
}

MKL_INT q_patch_grid_num_el(q_patch_t *q_patch) {
    return (q_patch->n_xi)*(q_patch->n_eta);
}

void q_patch_evaluate_M_p(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, rd_mat_t *x, rd_mat_t *y) {
    rd_mat_shape(x, xi.rows, xi.columns);
    rd_mat_shape(y, xi.rows, xi.columns);

    q_patch->M_p.M_p_handle(xi, eta, x, y, q_patch->M_p.extra_param);
}
void q_patch_evaluate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t *J_vals) {
    q_patch->J.J_handle(v, J_vals, q_patch->J.extra_param);
}

void q_patch_xi_mesh(q_patch_t *q_patch, rd_mat_t *xi_mesh_vals) {
    rd_mat_shape(xi_mesh_vals, q_patch->n_xi, 1);
    rd_linspace(q_patch->xi_start, q_patch->xi_end, q_patch->n_xi, xi_mesh_vals);
}

void q_patch_eta_mesh(q_patch_t *q_patch, rd_mat_t *eta_mesh_vals) {
    rd_mat_shape(eta_mesh_vals, q_patch->n_eta, 1);
    rd_linspace(q_patch->eta_start, q_patch->eta_end, q_patch->n_eta, eta_mesh_vals);
}

void q_patch_xi_eta_mesh(q_patch_t *q_patch, rd_mat_t *XI_vals, rd_mat_t *ETA_vals) {
    double xi_mesh_data[q_patch->n_xi];
    rd_mat_t xi_mesh = rd_mat_init_no_shape(xi_mesh_data);
    q_patch_xi_mesh(q_patch, &xi_mesh);

    double eta_mesh_data[q_patch->n_eta];
    rd_mat_t eta_mesh = rd_mat_init_no_shape(eta_mesh_data);
    q_patch_eta_mesh(q_patch, &eta_mesh);

    rd_meshgrid(xi_mesh, eta_mesh, XI_vals, ETA_vals);
}

void q_patch_convert_to_XY(q_patch_t *q_patch, rd_mat_t XI, rd_mat_t ETA, rd_mat_t *X_vals, rd_mat_t *Y_vals) {
    q_patch_evaluate_M_p(q_patch, XI, ETA, X_vals, Y_vals);
}

void q_patch_xy_mesh(q_patch_t *q_patch, rd_mat_t *X_vals, rd_mat_t *Y_vals) {
    double XI_data[q_patch_grid_num_el(q_patch)];
    double ETA_data[q_patch_grid_num_el(q_patch)];
    rd_mat_t XI = rd_mat_init_no_shape(XI_data);
    rd_mat_t ETA = rd_mat_init_no_shape(ETA_data);

    q_patch_xi_eta_mesh(q_patch, &XI, &ETA);
    q_patch_convert_to_XY(q_patch, XI, ETA, X_vals, Y_vals);
}

void q_patch_evaluate_f(q_patch_t *q_patch, scalar_func_2D_t f) {
    double X_data[q_patch_grid_num_el(q_patch)];
    double Y_data[q_patch_grid_num_el(q_patch)];
    rd_mat_t X = rd_mat_init_no_shape(X_data);
    rd_mat_t Y = rd_mat_init_no_shape(Y_data);

    q_patch_xy_mesh(q_patch, &X, &Y);
    for (int i = 0; i < q_patch->n_eta*q_patch->n_xi; i++) {
        q_patch->f_XY->mat_data[i] = f(X_data[i], Y_data[i]);
    }
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
            xi->mat_data[i] = q_patch->xi_start;
        }
        else if (fabs(xi->mat_data[i]-q_patch->xi_end) < q_patch->eps_xi_eta) {
            xi->mat_data[i] = q_patch->xi_end;
        }
        if (fabs(eta->mat_data[i]-q_patch->eta_start) < q_patch->eps_xi_eta) {
            eta->mat_data[i] = q_patch->eta_start;
        }
        else if (fabs(eta->mat_data[i]-q_patch->eta_end) < q_patch->eps_xi_eta) {
            eta->mat_data[i] = q_patch->eta_end;
        }
    }
}

inverse_M_p_return_type_t q_patch_inverse_M_p(q_patch_t *q_patch, double x, double y, rd_mat_t* initial_guesses_xi, rd_mat_t* initial_guesses_eta) {
    // global data for the case no initial guesses are given
    int N = 20;
    int N_segment = ceil(N/4);
    double initial_guesses_xi_data[N_segment*4];
    double initial_guesses_eta_data[N_segment*4];
    rd_mat_t initial_guesses_xi_mat;
    rd_mat_t initial_guesses_eta_mat;

    if (initial_guesses_xi == NULL || initial_guesses_eta == NULL) {        
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

    int converged;

    for (MKL_INT k = 0; k < initial_guesses_xi->rows; k++) {
        converged = 0;
        v.mat_data[0] = initial_guesses_xi->mat_data[k];
        v.mat_data[1] = initial_guesses_eta->mat_data[k];

        //newton solve with max iterations 100 and error tolerance given by q_patch
        for (MKL_INT i = 0; i < 100; i++) {
            // evaluates difference between solutions in real space
            q_patch_evaluate_M_p(q_patch, rd_mat_init(v.mat_data, 1, 1), rd_mat_init(v.mat_data+1, 1, 1), &M_p_v_x, &M_p_v_y);
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
            q_patch_evaluate_J(q_patch, v, &J_addr);
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

    shift_idx_mesh(&interpol_xi_j_mesh, 0, q_patch->n_xi-1);
    shift_idx_mesh(&interpol_eta_j_mesh, 0, q_patch->n_eta-1);

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
            interpol_xi_idxs[i] = sub2ind(q_patch->n_eta, q_patch->n_xi, (sub_t) {interpol_eta_j_mesh_data[horz_idx], interpol_xi_j_mesh_data[i]});
        }
        vdPackV(M, q_patch->f_XY->mat_data, interpol_xi_idxs, interpol_val_data);
        interpol_xi_exact_data[horz_idx] = barylag(interpol_xi_mesh, interpol_val, xi);
    }

    double f_xy = barylag(interpol_eta_mesh, interpol_xi_exact, eta);
    return (locally_compute_return_type_t) {f_xy, 1};
}

typedef struct w_param {
    double R;
    double theta_0;
} w_param_t;

void evaluate_w_1D(rd_mat_t theta, rd_mat_t *w_vals, w_1D_t w_1D, w_param_t w_param) {
    rd_mat_shape(w_vals, theta.rows, theta.columns);
    double theta_transform_data[theta.rows*theta.columns];
    vdSubI(theta.rows*theta.columns, theta.mat_data, 1, &(w_param.theta_0), 0, theta_transform_data, 1);
    vdDivI(theta.rows*theta.columns, theta_transform_data, 1, &(w_param.R), 0, theta_transform_data, 1);

    rd_mat_t theta_transform = rd_mat_init(theta_transform_data, theta.rows, theta.columns);
    w_1D(theta_transform, w_vals);
}

double compute_xi_corner(q_patch_t *main_patch, q_patch_t *window_patch, bool window_fix_xi, double window_fixed_edge, bool window_patch_right);
w_param_t apply_w_normalization_window(q_patch_t *main_patch, w_param_t main_w_param, q_patch_t *window_patch, double window_xi_corner, bool up_down);

void q_patch_apply_w_normalization_xi_right(q_patch_t *main_patch, q_patch_t *window_patch) {
    double main_xi_corner = compute_xi_corner(main_patch, window_patch, true, window_patch->xi_end, true);
    double window_xi_corner = compute_xi_corner(window_patch, main_patch, true, main_patch->xi_end, false);

    w_param_t main_w_param = {main_xi_corner-main_patch->xi_end, main_patch->xi_end};

    apply_w_normalization_window(main_patch, main_w_param, window_patch, window_xi_corner, false);

}

double compute_xi_corner(q_patch_t *main_patch, q_patch_t *window_patch, bool window_fix_xi, double window_fixed_edge, bool window_patch_right) {
    double window_xi_edge;
    double window_eta_edge;
    if(window_fix_xi) {
        window_xi_edge = window_fixed_edge;
        window_eta_edge = window_patch->eta_start;
    } else {
        window_xi_edge = window_patch->xi_start;
        window_eta_edge = window_fixed_edge;
    }

    double main_xi_corner;
    if(window_patch_right) {
        main_xi_corner = main_patch->xi_end;
    }else {
        main_xi_corner = main_patch->xi_start;
    }

    bool first_iter = true;

    //preallocation of data for for loop
    rd_mat_t window_xi_edge_mat = rd_mat_init(&window_xi_edge, 1, 1);
    rd_mat_t window_eta_edge_mat = rd_mat_init(&window_eta_edge, 1, 1);

    double window_x_edge;
    double window_y_edge;
    rd_mat_t window_x_edge_mat = rd_mat_init(&window_x_edge, 1, 1);
    rd_mat_t window_y_edge_mat = rd_mat_init(&window_y_edge, 1, 1);

    double main_xi_edge;
    double main_eta_edge;
    rd_mat_t main_xi_edge_mat = rd_mat_init(&main_xi_edge, 1, 1);
    rd_mat_t main_eta_edge_mat = rd_mat_init(&main_eta_edge, 1, 1);
    while (true) {
        q_patch_evaluate_M_p(window_patch, window_xi_edge_mat, window_eta_edge_mat, &window_x_edge_mat, &window_y_edge_mat);
        
        inverse_M_p_return_type_t main_xi_eta;
        if(first_iter) {
            main_xi_eta = q_patch_inverse_M_p(main_patch, window_x_edge, window_y_edge, NULL, NULL);
            first_iter = false;
        } else {
            main_xi_eta = q_patch_inverse_M_p(main_patch, window_x_edge, window_y_edge, &main_xi_edge_mat, &main_eta_edge_mat);
        }
        main_xi_edge = main_xi_eta.xi;
        main_eta_edge = main_xi_eta.eta;
        if(!main_xi_eta.converged) {
            printf("Nonconvergence in computing boundary mesh values!!!\n");
            break;
        }

        if ((main_xi_edge < main_xi_corner && window_patch_right) || (main_xi_edge > main_xi_corner && !window_patch_right)) {
            main_xi_corner = main_xi_edge;
        }
        if (main_eta_edge > main_patch->eta_end || main_eta_edge < main_patch->eta_start) {
            break;
        }
        if (window_fix_xi) {
            window_eta_edge += window_patch->h_eta;
        } else {
            window_xi_edge += window_patch->h_xi;
        }
    }

    return main_xi_corner;
}

MKL_INT xi_overlap_mesh_num_el(q_patch_t *main_patch, double xi_corner, bool window_patch_right) {
    if (window_patch_right) {
        MKL_INT xi_corner_j = ceil((xi_corner-main_patch->xi_start)/main_patch->h_xi);
        return (main_patch->n_xi-xi_corner_j) * (main_patch->n_eta);
    } else {
        MKL_INT xi_corner_j = floor((xi_corner-main_patch->xi_start)/main_patch->h_xi);
        return (xi_corner_j+1) * (main_patch->n_eta);
    }
}

void compute_xi_overlap_mesh(q_patch_t *main_patch, double xi_corner, bool window_patch_right, rd_mat_t *XI_overlap, rd_mat_t *ETA_overlap, ri_mat_t *XI_j, ri_mat_t *ETA_j) {
    if (window_patch_right) {
        MKL_INT xi_corner_j = ceil((xi_corner-main_patch->xi_start)/main_patch->h_xi);
        MKL_INT overlap_n_xi = main_patch->n_xi-xi_corner_j;

        MKL_INT xi_mesh_data[overlap_n_xi];
        ri_mat_t xi_mesh = ri_mat_init(xi_mesh_data, overlap_n_xi, 1);
        MKL_INT eta_mesh_data[main_patch->n_eta];
        ri_mat_t eta_mesh = ri_mat_init(eta_mesh_data, main_patch->n_eta, 1);

        ri_range(xi_corner_j, 1, main_patch->n_xi-1, &xi_mesh);
        ri_range(0, 1, main_patch->n_eta-1, &eta_mesh);

        ri_meshgrid(xi_mesh, eta_mesh, XI_j, ETA_j);

    }
    else {
        MKL_INT xi_corner_j = floor((xi_corner-main_patch->xi_start)/main_patch->h_xi);
        MKL_INT overlap_n_xi = xi_corner_j+1;

        MKL_INT xi_mesh_data[overlap_n_xi];
        ri_mat_t xi_mesh = ri_mat_init(xi_mesh_data, overlap_n_xi, 1);
        MKL_INT eta_mesh_data[main_patch->n_eta];
        ri_mat_t eta_mesh = ri_mat_init(eta_mesh_data, main_patch->n_eta, 1);

        ri_range(0, 1, xi_corner_j-1, &xi_mesh);
        ri_range(0, 1, main_patch->n_eta-1, &eta_mesh);

        ri_meshgrid(xi_mesh, eta_mesh, XI_j, ETA_j);
    }

    rd_mat_shape(XI_overlap, XI_j->rows, XI_j->columns);
    rd_mat_shape(ETA_overlap, XI_j->rows, XI_j->columns);
    for (int i = 0; i < XI_j->rows*XI_j->columns; i++) {
        XI_overlap->mat_data[i] = XI_j->mat_data[i] * main_patch->h_xi + main_patch->xi_start;
        ETA_overlap->mat_data[i] = ETA_j->mat_data[i] * main_patch->h_eta + main_patch->eta_start;
    }
}

void apply_w(q_patch_t *main_patch, rd_mat_t w_unnormalized, q_patch_t *window_patch, w_param_t window_w_param, rd_mat_t overlap_X, rd_mat_t overlap_Y, ri_mat_t overlap_XI_j, ri_mat_t overlap_ETA_j) {
    rd_mat_t *initial_guesses_xi = NULL;
    rd_mat_t *initial_guesses_eta = NULL;

    //preallocation for for loop
    double window_patch_xi;
    double window_patch_eta;
    
    rd_mat_t window_patch_xi_mat = rd_mat_init(&window_patch_xi, 1, 1);
    rd_mat_t window_patch_eta_mat = rd_mat_init(&window_patch_eta, 1, 1);

    double window_w;
    rd_mat_t window_w_mat = rd_mat_init(&window_w, 1, 1);

    for (int i = 0; i < overlap_X.rows; i++) {
        MKL_INT j_lst_data[overlap_X.columns];
        ri_mat_t j_lst = ri_mat_init(j_lst_data, overlap_X.columns, 1);
        if (i % 2 == 1) {
            ri_range(0, 1, overlap_X.columns-1, &j_lst);
        }else  {
            ri_range(overlap_X.columns-1, -1, 0, &j_lst);
        }

        for (int j_idx = 0; j_idx < overlap_X.columns; j_idx++) {
            int j = j_lst_data[j_idx];
            MKL_INT idx = sub2ind(overlap_X.rows, overlap_X.columns, (sub_t) {i, j});

            inverse_M_p_return_type_t window_patch_xi_eta = q_patch_inverse_M_p(window_patch, overlap_X.mat_data[idx], overlap_Y.mat_data[idx], initial_guesses_xi, initial_guesses_eta);
            window_patch_xi = window_patch_xi_eta.xi;
            window_patch_eta = window_patch_xi_eta.eta;

            if(window_patch_xi_eta.converged) {
                MKL_INT xi_j = overlap_XI_j.mat_data[idx];
                MKL_INT eta_j = overlap_ETA_j.mat_data[idx];

                MKL_INT f_idx = sub2ind(main_patch->n_eta, main_patch->n_xi, (sub_t) {eta_j, xi_j});

                evaluate_w_1D(window_patch_xi_mat, &window_w_mat, window_patch->w_1D, window_w_param);
                main_patch->f_XY->mat_data[f_idx] = main_patch->f_XY->mat_data[f_idx] * w_unnormalized.mat_data[idx] / (w_unnormalized.mat_data[idx] + window_w);
                initial_guesses_xi = &window_patch_xi_mat;
                initial_guesses_eta = &window_patch_eta_mat;
            } else {
                printf("Nonconvergence in computing C norm!!!");
            }
        }

    }
}

w_param_t apply_w_normalization_window(q_patch_t *main_patch, w_param_t main_w_param, q_patch_t *window_patch, double window_xi_corner, bool up_down) {
    w_param_t window_w_param;
    if(up_down) {
        window_w_param = (w_param_t) {window_xi_corner-window_patch->xi_start, window_patch->xi_start};
    } else {
        window_w_param = (w_param_t) {window_xi_corner-window_patch->xi_end, window_patch->xi_end};
    }

    MKL_INT window_xi_overlap_mesh_num_el = xi_overlap_mesh_num_el(window_patch, window_xi_corner, !up_down);

    MKL_INT XI_j_data[window_xi_overlap_mesh_num_el];
    MKL_INT ETA_j_data[window_xi_overlap_mesh_num_el];
    double XI_data[window_xi_overlap_mesh_num_el];
    double ETA_data[window_xi_overlap_mesh_num_el];
    ri_mat_t XI_j = ri_mat_init_no_shape(XI_j_data);
    ri_mat_t ETA_j = ri_mat_init_no_shape(ETA_j_data);
    rd_mat_t XI = rd_mat_init_no_shape(XI_data);
    rd_mat_t ETA = rd_mat_init_no_shape(ETA_data);

    compute_xi_overlap_mesh(window_patch, window_xi_corner, !up_down, &XI, &ETA, &XI_j, &ETA_j);
    
    double w_unnormalized_data[window_xi_overlap_mesh_num_el];
    rd_mat_t w_unnormalized = rd_mat_init_no_shape(w_unnormalized_data);
    evaluate_w_1D(XI, &w_unnormalized, window_patch->w_1D, window_w_param);

    double X_data[window_xi_overlap_mesh_num_el];
    rd_mat_t X = rd_mat_init_no_shape(X_data);
    double Y_data[window_xi_overlap_mesh_num_el];
    rd_mat_t Y = rd_mat_init_no_shape(Y_data);

    q_patch_convert_to_XY(window_patch, XI, ETA, &X, &Y);
    apply_w(window_patch, w_unnormalized, main_patch, main_w_param, X, Y, XI_j, ETA_j);

    return window_w_param;
}

void w_1D(rd_mat_t theta, rd_mat_t *w_1D_vals) {
    rd_mat_shape(w_1D_vals, theta.rows, theta.columns);
    MKL_INT size = theta.rows*theta.columns;
    double one = 1.0;
    cblas_dcopy(size, &one, 0, w_1D_vals->mat_data, 1);
    cblas_daxpy(size, -2, theta.mat_data, 1, w_1D_vals->mat_data, 1);
    cblas_dscal(size, 6, w_1D_vals->mat_data, 1);
    vdErfc(size, w_1D_vals->mat_data, w_1D_vals->mat_data);
    cblas_dscal(size, 0.5, w_1D_vals->mat_data, 1);
}