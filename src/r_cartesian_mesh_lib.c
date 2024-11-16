#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "r_cartesian_mesh_lib.h"
#include "num_linalg_lib.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double round_end_bound(double start_bound, double end_bound, double h) {
    return ceil((end_bound-start_bound)/h)*h + start_bound;
}

MKL_INT n_1D(double start_bound, double end_bound, double h) {
    return round((round_end_bound(start_bound, end_bound, h)-start_bound)/h) + 1;
}

MKL_INT n_total(double x_start, double x_end, double y_start, double y_end, double h) {
    return n_1D(x_start, x_end, h) * n_1D(y_start, y_end, h);
}

void inpolygon_mesh(rd_mat_t R_X, rd_mat_t R_Y, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_msk);

void r_cartesian_mesh_init(r_cartesian_mesh_obj_t *r_cartesian_mesh_obj, double x_start, double x_end, double y_start, double y_end, double h, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_interior, ri_mat_t *interior_idxs, rd_mat_t *f_XY) {
    r_cartesian_mesh_obj->x_start = x_start;
    r_cartesian_mesh_obj->y_start = y_start;

    r_cartesian_mesh_obj->x_end = round_end_bound(x_start, x_end, h);
    r_cartesian_mesh_obj->y_end = round_end_bound(y_start, y_end, h);

    r_cartesian_mesh_obj->h = h;
    
    r_cartesian_mesh_obj->n_x = n_1D(x_start, x_end, h);
    r_cartesian_mesh_obj->n_y = n_1D(y_start, y_end, h);

    double x_mesh_data[r_cartesian_mesh_obj->n_x];
    double y_mesh_data[r_cartesian_mesh_obj->n_y];
    rd_mat_t x_mesh = rd_mat_init(x_mesh_data, r_cartesian_mesh_obj->n_x, 1);
    rd_mat_t y_mesh = rd_mat_init(y_mesh_data, r_cartesian_mesh_obj->n_y, 1);
    rd_linspace(r_cartesian_mesh_obj->x_start, r_cartesian_mesh_obj->x_end, r_cartesian_mesh_obj->n_x, &x_mesh);
    rd_linspace(r_cartesian_mesh_obj->y_start, r_cartesian_mesh_obj->y_end, r_cartesian_mesh_obj->n_y, &y_mesh);

    double R_X_data[x_mesh.rows*y_mesh.rows];
    double R_Y_data[x_mesh.rows*y_mesh.rows];
    rd_mat_t R_X = rd_mat_init_no_shape(R_X_data);
    rd_mat_t R_Y = rd_mat_init_no_shape(R_Y_data);
    rd_meshgrid(x_mesh, y_mesh, &R_X, &R_Y);

    MKL_INT in_msk_data[r_cartesian_mesh_obj->n_x * r_cartesian_mesh_obj->n_y];
    ri_mat_t in_msk = ri_mat_init_no_shape(in_msk_data);
    inpolygon_mesh(R_X, R_Y, boundary_X, boundary_Y, &in_msk);
}

void inpolygon_mesh(rd_mat_t R_X, rd_mat_t R_Y, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_msk) {
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

    for(int row_idx = 0; row_idx < in_msk->rows; row_idx++) {
        bool in_interior = false;
        for (int col_idx = 0; col_idx < in_msk->columns; col_idx++) {
            MKL_INT idx = sub2ind(in_msk->rows, in_msk->columns, (sub_t) {row_idx, col_idx});
            if (in_msk->mat_data[idx] && !in_interior) {
                in_interior = true;
                in_msk->mat_data[idx] = 0;
            } else if (in_msk->mat_data[idx] && in_interior){
                in_interior = false;
            }
            else if (in_interior) {
                in_msk->mat_data[idx] = 1;
            }
        }
    }

    ri_print_matrix(*in_msk);
}