#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "r_cartesian_mesh_lib.h"
#include "num_linalg_lib.h"
#include "q_patch_lib.h"

double round_end_bound(double start_bound, double end_bound, double h) {
    return ceil((end_bound-start_bound)/h)*h + start_bound;
}

MKL_INT n_1D(double start_bound, double end_bound, double h) {
    return round((round_end_bound(start_bound, end_bound, h)-start_bound)/h) + 1;
}

MKL_INT r_cartesian_n_total(double x_start, double x_end, double y_start, double y_end, double h) {
    return n_1D(x_start, x_end, h) * n_1D(y_start, y_end, h);
}

MKL_INT inpolygon_mesh(rd_mat_t R_X, rd_mat_t R_Y, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_msk);

void r_cartesian_mesh_init(r_cartesian_mesh_obj_t *r_cartesian_mesh_obj, double x_start, double x_end, double y_start, double y_end, double h, rd_mat_t boundary_X, rd_mat_t boundary_Y, rd_mat_t *R_X, rd_mat_t *R_Y, ri_mat_t *in_interior, rd_mat_t *f_R) {
    r_cartesian_mesh_obj->x_start = x_start;
    r_cartesian_mesh_obj->y_start = y_start;

    r_cartesian_mesh_obj->x_end = round_end_bound(x_start, x_end, h);
    r_cartesian_mesh_obj->y_end = round_end_bound(y_start, y_end, h);

    r_cartesian_mesh_obj->h = h;
    
    r_cartesian_mesh_obj->n_x = n_1D(x_start, x_end, h);
    r_cartesian_mesh_obj->n_y = n_1D(y_start, y_end, h);

    r_cartesian_mesh_obj->R_X = R_X;
    r_cartesian_mesh_obj->R_Y = R_Y;
    r_cartesian_mesh_obj->in_interior = in_interior;

    rd_mat_shape(f_R, r_cartesian_mesh_obj->n_y, r_cartesian_mesh_obj->n_x);
    r_cartesian_mesh_obj->f_R = f_R;

    double x_mesh_data[r_cartesian_mesh_obj->n_x];
    double y_mesh_data[r_cartesian_mesh_obj->n_y];
    rd_mat_t x_mesh = rd_mat_init(x_mesh_data, r_cartesian_mesh_obj->n_x, 1);
    rd_mat_t y_mesh = rd_mat_init(y_mesh_data, r_cartesian_mesh_obj->n_y, 1);
    rd_linspace(r_cartesian_mesh_obj->x_start, r_cartesian_mesh_obj->x_end, r_cartesian_mesh_obj->n_x, &x_mesh);
    rd_linspace(r_cartesian_mesh_obj->y_start, r_cartesian_mesh_obj->y_end, r_cartesian_mesh_obj->n_y, &y_mesh);

    rd_meshgrid(x_mesh, y_mesh, R_X, R_Y);
    inpolygon_mesh(*R_X, *R_Y, boundary_X, boundary_Y, in_interior);
}

void r_cartesian_mesh_interpolate_patch(r_cartesian_mesh_obj_t *r_cartesian_mesh_obj, q_patch_t *q_patch, MKL_INT M) {
    double bound_X_data[q_patch_boundary_mesh_num_el(q_patch)];
    double bound_Y_data[q_patch_boundary_mesh_num_el(q_patch)];
    rd_mat_t bound_X = rd_mat_init_no_shape(bound_X_data);
    rd_mat_t bound_Y = rd_mat_init_no_shape(bound_Y_data);

    q_patch_boundary_mesh_xy(q_patch, false, &bound_X, &bound_Y);
    
    MKL_INT in_patch_data[r_cartesian_mesh_obj->n_x * r_cartesian_mesh_obj->n_y];
    ri_mat_t in_patch = ri_mat_init_no_shape(in_patch_data);
    MKL_INT n_in_patch = inpolygon_mesh(*(r_cartesian_mesh_obj->R_X), *(r_cartesian_mesh_obj->R_Y), bound_X, bound_Y, &in_patch);

    for (int i = 0; i < r_cartesian_mesh_obj->n_x * r_cartesian_mesh_obj->n_y; i++) {
        if (r_cartesian_mesh_obj->in_interior->mat_data[i] && in_patch_data[i]) {
            n_in_patch -= 1;
        }
        in_patch_data[i] = in_patch_data[i] && !(r_cartesian_mesh_obj->in_interior->mat_data[i]);        
    }

    MKL_INT r_patch_idxs_data[n_in_patch];
    ri_mat_t r_patch_idxs = ri_mat_init(r_patch_idxs_data, n_in_patch, 1);
    MKL_INT curr_idx = 0;
    for (int i = 0; i < r_cartesian_mesh_obj->n_x * r_cartesian_mesh_obj->n_y; i++) {
        if (in_patch.mat_data[i]) {
            r_patch_idxs.mat_data[curr_idx] = i;
            in_patch.mat_data[i] = curr_idx+1;
            curr_idx += 1;
        }
    }

    MKL_INT n_patch_grid = q_patch_grid_num_el(q_patch);
    double patch_XI_data[n_patch_grid];
    double patch_ETA_data[n_patch_grid];
    rd_mat_t patch_XI = rd_mat_init_no_shape(patch_XI_data);
    rd_mat_t patch_ETA = rd_mat_init_no_shape(patch_ETA_data);
    q_patch_xi_eta_mesh(q_patch, &patch_XI, &patch_ETA);

    double patch_X_data[n_patch_grid];
    double patch_Y_data[n_patch_grid];
    rd_mat_t patch_X = rd_mat_init_no_shape(patch_X_data);
    rd_mat_t patch_Y = rd_mat_init_no_shape(patch_Y_data);
    q_patch_convert_to_XY(q_patch, patch_XI, patch_ETA, &patch_X, &patch_Y);

    MKL_INT floor_X_j[n_patch_grid];
    MKL_INT ceil_X_j[n_patch_grid];
    MKL_INT floor_Y_j[n_patch_grid];
    MKL_INT ceil_Y_j[n_patch_grid];

    for (int i = 0; i < n_patch_grid; i++) {
        floor_X_j[i] = (MKL_INT) floor((patch_X.mat_data[i] - r_cartesian_mesh_obj->x_start)/r_cartesian_mesh_obj->h);
        ceil_X_j[i] = (MKL_INT) ceil((patch_X.mat_data[i] - r_cartesian_mesh_obj->x_start)/r_cartesian_mesh_obj->h);
        floor_Y_j[i] = (MKL_INT) floor((patch_Y.mat_data[i] - r_cartesian_mesh_obj->y_start)/r_cartesian_mesh_obj->h);
        ceil_Y_j[i] = (MKL_INT) ceil((patch_Y.mat_data[i] - r_cartesian_mesh_obj->y_start)/r_cartesian_mesh_obj->h);
    }

    double P_xi[n_in_patch];
    double P_eta[n_in_patch];

    for (int i = 0; i < n_in_patch; i++) {
        P_xi[i] = NAN;
        P_eta[i] = NAN;
    }

    for (int i = 0; i < n_patch_grid; i++) {
        double neighbors_X[4] = {floor_X_j[i], floor_X_j[i], ceil_X_j[i], ceil_X_j[i]};
        double neighbors_Y[4] = {floor_Y_j[i], ceil_Y_j[i], floor_Y_j[i], ceil_Y_j[i]};

        for (int j = 0; j < 4; j++) {
            double neighbor_X = neighbors_X[j];
            double neighbor_Y = neighbors_Y[j];

            if (neighbor_X > r_cartesian_mesh_obj->n_x-1 || neighbor_X < 0 || neighbor_Y > r_cartesian_mesh_obj->n_y-1 || neighbor_Y < 0) {
                continue;
            }
            MKL_INT patch_idx = sub2ind(r_cartesian_mesh_obj->n_y, r_cartesian_mesh_obj->n_x, (sub_t) {neighbor_Y, neighbor_X});
            
            if(in_patch.mat_data[patch_idx] != 0 && isnan(P_xi[in_patch.mat_data[patch_idx]-1])) {
                double neighbor_X_coord = neighbor_X*r_cartesian_mesh_obj->h+r_cartesian_mesh_obj->x_start;
                double neighbor_Y_coord = neighbor_Y*r_cartesian_mesh_obj->h+r_cartesian_mesh_obj->y_start;
                rd_mat_t initial_guess_xi = rd_mat_init(patch_XI_data+i, 1, 1);
                rd_mat_t initial_guess_eta = rd_mat_init(patch_ETA_data+i, 1, 1);
                inverse_M_p_return_type_t xi_eta = q_patch_inverse_M_p(q_patch, neighbor_X_coord, neighbor_Y_coord, &initial_guess_xi, &initial_guess_eta);

                if (xi_eta.converged) {
                    P_xi[in_patch.mat_data[patch_idx]-1] = xi_eta.xi;
                    P_eta[in_patch.mat_data[patch_idx]-1] = xi_eta.eta;
                }
                else {
                    printf("Nonconvergence in interpolation!!!");
                }
            }
        }
    }

    MKL_INT nan_count;
    bool first_iter = true;
    while (true) {
        if (!first_iter && nan_count == 0) {
            break;
        }
        nan_count = 0;
        for (int i = 0; i < n_in_patch; i++) {
            if (isnan(P_xi[i])) {
                bool is_touched = false;
                sub_t idx = ind2sub(r_cartesian_mesh_obj->n_y, r_cartesian_mesh_obj->n_x, r_patch_idxs_data[i]);

                MKL_INT neighbor_shift_x[8] = {1,-1, 0, 0, 1, 1, -1, -1};
                MKL_INT neighbor_shift_y[8] = {0, 0, -1, 1, 1, -1, 1, -1};
                for (int j = 0; j < 8; j++) {
                    MKL_INT neighbor_i = idx.i + neighbor_shift_y[j];
                    MKL_INT neighbor_j = idx.j + neighbor_shift_x[j];
                    
                    if (neighbor_j > r_cartesian_mesh_obj->n_x-1 || neighbor_j < 0 || neighbor_i > r_cartesian_mesh_obj->n_y-1 || neighbor_i < 0) {
                        continue;
                    }

                    MKL_INT neighbor = sub2ind(r_cartesian_mesh_obj->n_y, r_cartesian_mesh_obj->n_x, (sub_t) {neighbor_i, neighbor_j});
                    if (in_patch.mat_data[neighbor] != 0 && !isnan(P_xi[in_patch.mat_data[neighbor]-1])) {
                        double idx_x_coord = idx.j*r_cartesian_mesh_obj->h + r_cartesian_mesh_obj->x_start;
                        double idx_y_coord = idx.i*r_cartesian_mesh_obj->h + r_cartesian_mesh_obj->y_start;
                        rd_mat_t initial_guess_xi = rd_mat_init(P_xi + in_patch.mat_data[neighbor] - 1, 1, 1);
                        rd_mat_t initial_guess_eta = rd_mat_init(P_eta + in_patch.mat_data[neighbor] - 1, 1, 1);
                        inverse_M_p_return_type_t xi_eta = q_patch_inverse_M_p(q_patch, idx_x_coord, idx_y_coord, &initial_guess_xi, &initial_guess_eta);

                        if (xi_eta.converged) {
                            P_xi[i] = xi_eta.xi;
                            P_eta[i] = xi_eta.eta;
                            is_touched = true;
                        }
                        else {
                            printf("Nonconvergence in interpolation!!!");
                        }
                    }
                }
                if (!is_touched) {
                    nan_count += 1;
                }
            }
        }
        first_iter = false;
    }

    double f_R_patch[n_in_patch];
    for (int i = 0; i < n_in_patch; i++) {
        double xi_point = P_xi[in_patch.mat_data[r_patch_idxs.mat_data[i]]-1];
        double eta_point = P_eta[in_patch.mat_data[r_patch_idxs.mat_data[i]]-1];

        locally_compute_return_type_t f_locally_compute = q_patch_locally_compute(q_patch, xi_point, eta_point, M);
        if(f_locally_compute.in_range) {
            f_R_patch[i] = f_locally_compute.f_xy;
        }
    }
    
    for (int i = 0; i < n_in_patch; i++) {
        r_cartesian_mesh_obj->f_R->mat_data[r_patch_idxs.mat_data[i]] += f_R_patch[i];
    }
}

MKL_INT inpolygon_mesh(rd_mat_t R_X, rd_mat_t R_Y, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_msk) {
    ri_mat_shape(in_msk, R_X.rows, R_X.columns);
    memset(in_msk->mat_data, 0, in_msk->rows*in_msk->columns*sizeof(MKL_INT));

    MKL_INT n_edges = boundary_X.rows-1;

    rd_mat_t boundary_x_edge_1 = rd_mat_init(boundary_X.mat_data, n_edges, 1);
    rd_mat_t boundary_x_edge_2 = rd_mat_init(boundary_X.mat_data+1, n_edges, 1);
    rd_mat_t boundary_y_edge_1 = rd_mat_init(boundary_Y.mat_data, n_edges, 1);
    rd_mat_t boundary_y_edge_2 = rd_mat_init(boundary_Y.mat_data+1, n_edges, 1);

    MKL_INT boundary_idx_data[n_edges];
    ri_mat_t boundary_idx = ri_mat_init(boundary_idx_data, n_edges, 1);
    ri_range(0, 1, n_edges-1, &boundary_idx);

    double x_start = R_X.mat_data[0];
    double y_start = R_Y.mat_data[0];
    double h_x = R_X.mat_data[R_X.rows] - R_X.mat_data[0];
    double h_y = R_Y.mat_data[1] - R_Y.mat_data[0];

    double boundary_y_j_data[boundary_Y.rows];
    vdSubI(boundary_Y.rows, boundary_Y.mat_data, 1, &y_start, 0, boundary_y_j_data, 1);
    vdDivI(boundary_Y.rows, boundary_y_j_data, 1, &h_y, 0, boundary_y_j_data, 1);

    rd_mat_t boundary_y_edge_1_j = rd_mat_init(boundary_y_j_data, n_edges, 1);
    rd_mat_t boundary_y_edge_2_j = rd_mat_init(boundary_y_j_data+1, n_edges, 1);

    MKL_INT intersection_edges_msk[n_edges];
    MKL_INT n_intersection_edges = 0;
    for (int i = 0; i < n_edges; i++ ) {
        if(fabs(boundary_y_edge_1_j.mat_data[i] - round(boundary_y_edge_1_j.mat_data[i])) < __DBL_EPSILON__) {
            boundary_y_edge_1_j.mat_data[i] = round(boundary_y_edge_1_j.mat_data[i]);
        }
        if(fabs(boundary_y_edge_2_j.mat_data[i] - round(boundary_y_edge_2_j.mat_data[i])) < __DBL_EPSILON__) {
            boundary_y_edge_2_j.mat_data[i] = round(boundary_y_edge_2_j.mat_data[i]);
        }

        intersection_edges_msk[i] = floor(boundary_y_edge_1_j.mat_data[i]) != floor(boundary_y_edge_2_j.mat_data[i]);
        if (floor(boundary_y_edge_1_j.mat_data[i]) != floor(boundary_y_edge_2_j.mat_data[i])) {
            n_intersection_edges += 1;
        }
    }
    MKL_INT intersection_idxs_data[n_intersection_edges];
    MKL_INT curr_idx = 0;
    for (int i = 0; i < n_edges; i++) {
        if(intersection_edges_msk[i]) {
            intersection_idxs_data[curr_idx] = boundary_idx_data[i];
            curr_idx += 1;
        }
    }
    
    for (int i = 0; i < n_intersection_edges; i++) {
        MKL_INT intersection_idx = intersection_idxs_data[i];

        double x_edge_1 = boundary_x_edge_1.mat_data[intersection_idx];
        double x_edge_2 = boundary_x_edge_2.mat_data[intersection_idx];
        double y_edge_1 = boundary_y_edge_1.mat_data[intersection_idx];
        double y_edge_2 = boundary_y_edge_2.mat_data[intersection_idx];
        MKL_INT y_edge_1_j = (MKL_INT) floor(boundary_y_edge_1_j.mat_data[intersection_idx]);
        MKL_INT y_edge_2_j = (MKL_INT) floor(boundary_y_edge_2_j.mat_data[intersection_idx]);

        MKL_INT intersection_mesh_length = MAX(y_edge_1_j, y_edge_2_j) - MIN(y_edge_1_j, y_edge_2_j);
        MKL_INT intersection_mesh_y_j_data[intersection_mesh_length];
        ri_mat_t intersection_mesh_y_j = ri_mat_init(intersection_mesh_y_j_data, intersection_mesh_length, 1);
        ri_range(MIN(y_edge_1_j, y_edge_2_j)+1, 1, MAX(y_edge_1_j, y_edge_2_j), &intersection_mesh_y_j);

        for (int j = 0; j < intersection_mesh_length; j++) {
            double intersection_y = intersection_mesh_y_j_data[j] * h_y + y_start;
            double intersection_x = x_edge_1 + (x_edge_2-x_edge_1)*(intersection_y-y_edge_1)/(y_edge_2-y_edge_1);

            MKL_INT mesh_intersection_idx = sub2ind(in_msk->rows, in_msk->columns, (sub_t) {(MKL_INT) round((intersection_y-y_start)/h_y), floor((intersection_x-x_start)/h_x)});
            in_msk->mat_data[mesh_intersection_idx] = !in_msk->mat_data[mesh_intersection_idx];
        }
    }

    MKL_INT n_points_interior = 0;
    for(int row_idx = 0; row_idx < in_msk->rows; row_idx++) {
        bool in_interior = false;
        for (int col_idx = 0; col_idx < in_msk->columns; col_idx++) {
            MKL_INT idx = sub2ind(in_msk->rows, in_msk->columns, (sub_t) {row_idx, col_idx});
            if (in_msk->mat_data[idx] && !in_interior) {
                in_interior = true;
                in_msk->mat_data[idx] = 0;
            } else if (in_msk->mat_data[idx] && in_interior){
                in_interior = false;
                n_points_interior += 1;
            }
            else if (in_interior) {
                in_msk->mat_data[idx] = 1;
                n_points_interior += 1;
            }
        }
    }

    return n_points_interior;
}